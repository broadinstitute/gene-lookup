import base64
import datetime
from dotenv import load_dotenv
import json
import pandas as pd
import requests
import os
import urllib.request
from tqdm import tqdm
import re

from annotation_utils.cache_utils import cache_data_table

load_dotenv()

OUTPUT_COLUMNS = [
    'chrom',
    'start',
    'end',
    'mim_number',
    'phenotype_map_method',
    'phenotype_mim_number',
    'phenotype_inheritance',
    'gene_symbols',
    'gene_id',
    'gene_description',
    'phenotype_description',
    'mouse_gene_id',
    'comments',
]

MAX_GENE_SIZE = 5*10**6  # 5 Mbases (dystrophin is 2.3Mb)

"""
Example row:
$1                      chrom : 1
$2                      start : 7784284
$3                        end : 7845180
$4                 mim_number : 603427
$5       phenotype_mim_number : 616882
$6   phenotypic_series_number :
$7      phenotype_inheritance : Autosomal dominant
$8               gene_symbols : FASPS3, PER3
$9                    gene_id : ENSG00000049246
$10          gene_description : Period circadian regulator 3
$11     phenotype_description : ?Advanced sleep phase syndrome, familial, 3
$12              date_created :
$13              date_updated :
$14             mouse_gene_id : Per3 (MGI:1277134)
$15              oe_lof_upper : 8.4100e-01
$16                       pLI : 1.3661e-17
$17                     mis_z : -1.0855e-01
$18                      text :
$19                  comments : mutation identified in 1 FASPS3 family
"""

OMIM_PHENOTYPE_MAP_METHOD_CHOICES = {
    1: 'the disorder is placed on the map based on its association with a gene, but the underlying defect is not known.',
    2: 'the disorder has been placed on the map by linkage; no mutation has been found.',
    3: 'the molecular basis for the disorder is known; a mutation has been found in the gene.',
    4: 'a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype.',
}


def download_file(url, to_dir=".", force_download=False, verbose=True):
    """Download the given file and returns its local path.
     Args:
        url (string): HTTP or FTP url
     Returns:
        string: local file path
    """

    if not (url and url.startswith(("http://", "https://", "ftp://"))):
        raise ValueError("Invalid url: {}".format(url))

    local_file_path = os.path.join(to_dir, os.path.basename(url))
    if not force_download and os.path.isfile(local_file_path) and os.path.getsize(local_file_path) > 100:
        return local_file_path

    print(f"Downloading {url} to {local_file_path}")

    input_iter = urllib.request.urlopen(url)
    if verbose:
        input_iter = tqdm(input_iter, unit=" data" if url.endswith("gz") else " lines")

    with open(local_file_path, 'w') as f:
        f.writelines((line.decode('utf-8') for line in input_iter))

    input_iter.close()

    return local_file_path


def get_file_header(f):
    header_fields = None
    for i, line in enumerate(f):
        line = line.rstrip('\r\n')
        if line.startswith("# Chrom") and header_fields is None:
            header_fields = [c.lower().replace(' ', '_') for c in line.lstrip('# ').split('\t')]

            break
        elif not line or line.startswith("#"):
            continue
        elif line.startswith('This account is inactive') or line.startswith('This account has expired'):
            raise Exception(line)
        elif header_fields is None:
            raise ValueError("Header row not found in genemap2 file before line {}: {}".format(i, line))

    renamed_header_fields = []
    for header_field in header_fields:
        if header_field == "gene/locus_and_other_related_symbols":
            header_field = "gene_symbols"
            print(f"Renamed header field to gene_symbols")
        renamed_header_fields.append(header_field)

    return renamed_header_fields


def parse_genemap2_records(omim_line_fields):
    # skip commented rows
    if len(omim_line_fields) == 1:
        yield None

    else:
        # rename some of the fields
        output_record = {}
        output_record['chrom'] = omim_line_fields['chromosome'].replace("chr", "")
        output_record['start'] = int(omim_line_fields['genomic_position_start']) or 1  # 'or 1' replaces pos=0 with 1
        output_record['end'] = int(omim_line_fields['genomic_position_end']) or 1
        output_record['cyto'] = omim_line_fields['cyto_location']
        output_record['gene_id'] = omim_line_fields['ensembl_gene_id']
        output_record['mim_number'] = int(omim_line_fields['mim_number'])
        output_record['gene_symbols'] = ", ".join(sorted(set([s for s in [omim_line_fields.get('approved_symbol', '').strip()] + list(omim_line_fields['gene_symbols'].split(",")) if s])))
        output_record['gene_description'] = omim_line_fields['gene_name']
        output_record['comments'] = omim_line_fields['comments']
        output_record['mouse_gene_id'] = omim_line_fields['mouse_gene_symbol/id']

        phenotype_field = omim_line_fields['phenotypes'].strip()

        record_with_phenotype = None
        for phenotype_match in re.finditer("[\[{ ]*(.+?)[ }\]]*(, (\d{4,}))? \(([1-4])\)(, ([^;]+))?;?", phenotype_field):
            # Phenotypes example: "Langer mesomelic dysplasia, 249700 (3), Autosomal recessive; Leri-Weill dyschondrosteosis, 127300 (3), Autosomal dominant"

            record_with_phenotype = dict(output_record)  # copy
            record_with_phenotype["phenotype_description"] = phenotype_match.group(1)
            record_with_phenotype["phenotype_mim_number"] = int(phenotype_match.group(3)) if phenotype_match.group(3) else None
            record_with_phenotype["phenotype_map_method"] = int(phenotype_match.group(4))
            record_with_phenotype["phenotype_inheritance"] = phenotype_match.group(6) or None

            # basic checks
            if len(record_with_phenotype["phenotype_description"].strip()) == 0:
                raise ValueError("Empty phenotype description: {}".format(json.dumps(omim_line_fields)))

            if int(record_with_phenotype["phenotype_map_method"]) not in OMIM_PHENOTYPE_MAP_METHOD_CHOICES:
                raise ValueError("Unexpected value (%s) for phenotype_map_method: %s" % (
                    record_with_phenotype["phenotype_map_method"], phenotype_field))

            phenotype_map_method = OMIM_PHENOTYPE_MAP_METHOD_CHOICES[int(record_with_phenotype["phenotype_map_method"])]
            yield record_with_phenotype

        if record_with_phenotype is None:
            if len(phenotype_field) > 0:
                raise ValueError("No phenotypes found: {}".format(json.dumps(omim_line_fields)))
            else:
                yield output_record


@cache_data_table
def get_omim_table():
    """Retrieves the latest OMIM table"""

    omim_key = os.getenv("OMIM_KEY")
    if not omim_key:
        raise Exception("OMIM_KEY is not set in .env")

    file_path = download_file(f"https://data.omim.org/downloads/{omim_key}/genemap2.txt", force_download=True)
    #file_path = os.path.expanduser("~/code/omim-search-private/genemap2.txt")
    print(f'Parsing {file_path}')
    records = []
    with open(file_path) as f:
        header_fields = get_file_header(f)

        for line in tqdm(f, unit=" records"):
            omim_line_fields = dict(zip(header_fields, line.rstrip('\r\n').split('\t')))
            for record in parse_genemap2_records(omim_line_fields):
                if record is None:
                    continue

                records.append(record)

    omim_df = pd.DataFrame(records)

    print(f"Got {len(omim_df[omim_df['gene_id'].notna() & omim_df['phenotype_mim_number'].notna()].gene_id.unique()):,d} unique genes from OMIM")

    #print(omim_df.columns)
    #print(omim_df["locus"])
    #omim_df = omim_df[omim_df["locus_size"] < MAX_GENE_SIZE]
    #omim_df = omim_df[omim_df["start"] > 1]

    """
    Autosomal recessive                                          3333
    Autosomal dominant                                           2480
    X-linked recessive                                            217
    Autosomal dominant, Autosomal recessive                       209
    X-linked                                                       75
    X-linked dominant                                              69
    Somatic mutation, Autosomal dominant                           60
    Multifactorial                                                 31
    Multifactorial, Autosomal dominant, Autosomal recessive        12
    Digenic recessive, Autosomal recessive                         11
    Somatic mutation                                               10
    Autosomal dominant, Multifactorial                              9
    Digenic recessive                                               8
    Isolated cases, Autosomal dominant                              7
    Digenic dominant, Autosomal recessive                           6
    X-linked dominant, X-linked recessive                           5
    Isolated cases                                                  5
    Digenic dominant                                                5
    ?Autosomal dominant                                             4
    Digenic dominant, Autosomal dominant, Autosomal recessive       4
    Somatic mutation, Autosomal dominant, Autosomal recessive       4
    Digenic dominant, Autosomal dominant                            3
    X-linked, Isolated cases, Multifactorial                        3
    Y-linked                                                        3
    Multifactorial, Autosomal recessive                             2
    Somatic mosaicism, Autosomal recessive                          2
    Pseudoautosomal dominant                                        2
    Pseudoautosomal recessive                                       2
    Isolated cases, Somatic mutation                                1
    Mitochondrial                                                   1
    ?Autosomal dominant, Autosomal recessive                        1
    X-linked, Somatic mosaicism                                     1
    Somatic mosaicism, Autosomal dominant                           1
    X-linked dominant, Somatic mosaicism                            1
    """    

    # rename inheritance
    omim_df["phenotype_inheritance"] = omim_df["phenotype_inheritance"].map({
        "Autosomal recessive": "AR",
        "Autosomal dominant": "AD",
        "X-linked recessive": "XR",
        "X-linked dominant": "XD",
        "X-linked": "XR",
        "Autosomal dominant, Autosomal recessive": "AD/AR",
        "X-linked dominant, X-linked recessive": "XD/XR",
        "Somatic mutation, Autosomal dominant": "Somatic/AD",
        "Somatic mutation, Autosomal recessive": "Somatic/AR",
        "Somatic mutation, Autosomal dominant, Autosomal recessive": "Somatic/AD/AR",
        "Somatic mosaicism, Autosomal dominant": "Somatic/AD",
        "Somatic mosaicism, Autosomal recessive": "Somatic/AR",
        "Somatic mutation": "Somatic",
        "Mitochondrial": "MITO",
    })

    return omim_df[OUTPUT_COLUMNS]


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_omim_table()
    print(df)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_path = f"omim_{timestamp}.tsv"
    df.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(df)} records to {output_path}")


    print(df["phenotype_inheritance"].value_counts())
"""
5/19/2025 OMIM Table:

Columns:

0: chrom
1: start
2: end
3: mim_number
4: phenotype_mim_number
5: phenotypic_series_number
6: phenotype_inheritance
7: gene_symbols
8: gene_id
9: gene_description
10: phenotype_description
11: date_created
12: date_updated
13: mouse_gene_id
14: oe_lof_upper
15: pLI
16: mis_z
17: text
18: comments


(19525, 19)

Example:

chrom                                                   1
start                                             1013497
end                                               1014540
mim_number                                         147571
phenotype_mim_number                               616126
phenotypic_series_number
phenotype_inheritance                 Autosomal recessive
gene_symbols                  G1P2,  IFI15,  IMD38, ISG15
gene_id                                   ENSG00000187608
gene_description            ISG15 ubiquitin-like modifier
phenotype_description                 Immunodeficiency 38
date_created
date_updated
mouse_gene_id                         Isg15 (MGI:1855694)
oe_lof_upper                                        1.691
pLI                                               0.40527
mis_z                                           -0.044129
text
comments
"""
