import pandas as pd
import requests
import time

from annotation_utils.cache_utils import cache_data_table

URL_UK_PANEL_APP = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"
URL_AUSTRALIA_PANEL_APP = "https://panelapp-aus.org/api/v1/genes/"

@cache_data_table
def get_panel_app_table():
    """Download the PANEL app table from https://www.genenames.org/download/custom/ and return it as a pandas DataFrame"""
    rows = []
    for source_label, url in [("PanelApp UK", URL_UK_PANEL_APP), ("PanelApp Australia", URL_AUSTRALIA_PANEL_APP)]:
        data = {
            'next': url
        }

        page_i = 0
        while data.get('next'):  # page through the results
            page_i += 1
            url = data['next']

            time.sleep(2)  # sleep for 2 seconds to avoid being rate limited

            print(f"Retrieving page {page_i} of {url}")
            r = requests.get(url)
            if not r.ok:
                raise Exception(f"Failed to download {url}: {r}")

            data = r.json()

            for r in data["results"]:
                ensembl_genes = r["gene_data"]["ensembl_genes"]

                # gene id may not be specified for some results
                gene_id = ""
                if ensembl_genes and 'GRch38' in ensembl_genes:
                    gene_id = next(iter(ensembl_genes['GRch38'].values())).get("ensembl_id")

                rows.append({
                    "source": source_label,
                    "hgnc": r["gene_data"]["hgnc_id"],
                    "gene_name": r["gene_data"]["gene_name"],
                    "biotype": r["gene_data"]["biotype"],
                    "gene_id": gene_id,
                    "confidence": r["confidence_level"],
                    "penetrance": r["penetrance"],
                    "mode_of_pathogenicity": r["mode_of_pathogenicity"],
                    "mode_of_inheritance": r["mode_of_inheritance"],
                    "publications": ", ".join(r["publications"]).replace("\t", "  "),
                    "evidence": ", ".join(r["evidence"]).replace("\t", "  "),
                    "phenotypes": ", ".join(r["phenotypes"]).replace("\t", "  "),
                    "panel_name": r["panel"]["name"],
                })

        print(f"Retrieved {url}  total rows: {len(rows)}")

    df = pd.DataFrame(rows)

    # rename mode_of_inheritance
    """
    BIALLELIC, autosomal or pseudoautosomal                                                                                                    35467
    MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted                                                                                   11194
    MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown                                                                         6076
    BOTH monoallelic and biallelic, autosomal or pseudoautosomal                                                                                3571
    Unknown                                                                                                                                     3351
    X-LINKED: hemizygous mutation in males, biallelic mutations in females                                                                      2307
    X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)     1364
    BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal                      607
    MITOCHONDRIAL                                                                                                                                374
    Other                                                                                                                                        284
    MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)                                                  201
    MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed)                                                  138
    Other - please specifiy in evaluation comments                                                                                                34
    Other - please specify in evaluation comments                                                                                                 14
    """
    df["mode_of_inheritance"] = df["mode_of_inheritance"].map({
        "BIALLELIC, autosomal or pseudoautosomal": "AR",
        "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted": "AD",
        "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown": "AD",
        "BOTH monoallelic and biallelic, autosomal or pseudoautosomal": "AD/AR",
        "Unknown": "Unknown",
        "X-LINKED: hemizygous mutation in males, biallelic mutations in females": "XR",
        "X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)": "XD",
        "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal": "AD/AR (AR is more severe)",
        "MITOCHONDRIAL": "MITO",
        "Other": "Other",
        "MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)": "AD (maternally imprinted, paternal allele expressed)",
        "MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed)": "AD (paternally imprinted, maternal allele expressed)",
        "Other - please specifiy in evaluation comments": "Other",
        "Other - please specify in evaluation comments": "Other",
    })

    return df


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_panel_app_table()
    print(df)

    print(df["mode_of_inheritance"].value_counts())


"""
5/19/2025 PANEL APP Table:

Columns:

0: source
1: hgnc
2: gene_name
3: biotype
4: gene_id
5: confidence
6: penetrance
7: mode_of_pathogenicity
8: mode_of_inheritance
9: publications
10: evidence
11: phenotypes
12: panel_name

(200, 13)

Example:

source                                               PanelApp UK
hgnc                                                  HGNC:24641
gene_name                    chromosome 16 open reading frame 62
biotype                                           protein_coding
gene_id                                          ENSG00000103544
confidence                                                     2
penetrance                                                  None
mode_of_pathogenicity                                       None
mode_of_inheritance      BIALLELIC, autosomal or pseudoautosomal
publications                                            31712251
evidence                         Expert Review Amber, Literature
phenotypes                    3C/Ritscher-Schinzel-like syndrome
panel_name                             Chondrodysplasia punctata



total rows: 68957


source
PanelApp UK           35274
PanelApp Australia    33683


biotype
protein_coding                        68075
Mt_tRNA                                 161
snRNA                                    56
lincRNA                                  38
processed_transcript                     37
miRNA                                    33
snoRNA                                   28
antisense_RNA                            19
Mt_rRNA                                  18
IG_C_gene                                18
transcribed_unprocessed_pseudogene       10
processed_pseudogene                      4
unprocessed_pseudogene                    3
transcribed_unitary_pseudogene            2
vaultRNA                                  2
IG_V_gene                                 2
polymorphic_pseudogene                    1


confidence
3    46422
1    15401
2     6412
0      722


penetrance
Complete      12242
unknown         797
Incomplete      219


mode_of_pathogenicity
                                                                                                                                 50983
Other                                                                                                                             1454
Loss-of-function variants (as defined in pop up message) DO NOT cause this phenotype - please provide details in the comments      754
Other - please provide details in the comments                                                                                     358
Other - please provide details in the comments                                                                                      44
Loss-of-function variants (as defined in pop up message) DO NOT cause this phenotype -please provide details in the comments        31
                                                                                                                                    15
other - please provide details in the comments                                                                                       4
gain-of-function                                                                                                                     1
loss-of-function (truncating variants and curated list of variants)                                                                  1
FALSE                                                                                                                                1



mode_of_inheritance
BIALLELIC, autosomal or pseudoautosomal                                                                                                    35405
MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted                                                                                   11111
MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown                                                                         6100
                                                                                                                                            4125
BOTH monoallelic and biallelic, autosomal or pseudoautosomal                                                                                3587
Unknown                                                                                                                                     3358
X-LINKED: hemizygous mutation in males, biallelic mutations in females                                                                      2313
X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)     1370
BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal                      595
MITOCHONDRIAL                                                                                                                                312
Other                                                                                                                                        297
MONOALLELIC, autosomal or pseudoautosomal, maternally imprinted (paternal allele expressed)                                                  200
MONOALLELIC, autosomal or pseudoautosomal, paternally imprinted (maternal allele expressed)                                                  138
Other - please specifiy in evaluation comments                                                                                                32
Other - please specify in evaluation comments                                                                                                 14
                                                                                                                                           14


"""