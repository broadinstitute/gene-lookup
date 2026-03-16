import pandas as pd
import requests
from io import StringIO

from annotation_utils.cache_utils import cache_data_table
from annotation_utils.get_hgnc_table import get_hgnc_to_ensg_id_map

def _get_clingen_table(url):
    """Download one of the ClinGen .csv tables and return it as a pandas DataFrame

    Args:
        url (str): for example "https://search.clinicalgenome.org/kb/gene-validity/download"
    Return:
        pandas DataFrame
    """
    r = requests.get(url)
    if not r.ok:
        raise Exception(f"Failed to download {url}: {r}")

    table_contents = r.content.decode('UTF-8')
    lines = table_contents.split("\n")
    # Find the header line by looking for "GENE SYMBOL" rather than assuming a fixed position
    header_idx = None
    for i, line in enumerate(lines):
        if '"GENE SYMBOL"' in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Could not find header line containing 'GENE SYMBOL' in ClinGen download from {url}")
    # Skip separator line (all +++) after the header
    data_start = header_idx + 1
    if data_start < len(lines) and lines[data_start].startswith('"+++'):
        data_start += 1
    table_contents = "\n".join([lines[header_idx]] + lines[data_start:])
    return pd.read_csv(StringIO(table_contents))


@cache_data_table
def get_clingen_gene_disease_validity_table():
    "Download ClinGen gene-disease validity table and return it as a pandas DataFrame"
    return _get_clingen_table("https://search.clinicalgenome.org/kb/gene-validity/download")

@cache_data_table
def get_clingen_haploinsufficient_genes_table():
    """ClinGen’s Haploinsufficiency (HI) score ranges from 0 to 3, indicating how likely a gene is to cause disease due to loss-of-function (haploinsufficiency):

        0: No evidence for haploinsufficiency.
        1: Little evidence; haploinsufficiency is unlikely.
        2: Some evidence; haploinsufficiency is possible but not definitive.
        3: Sufficient evidence; haploinsufficiency is a known mechanism of disease.

    
    Values 30 and 40 in ClinGen's HI score system are historical and indicate haploinsufficient regions (e.g. microdeletion syndromes), rather than single genes:
        30: Some evidence a region is haploinsufficient.
        40: Strong evidence a region is haploinsufficient.
    """
    r = requests.get("https://search.clinicalgenome.org/kb/gene-dosage/download")
    if not r.ok:
        raise Exception(f"Failed to download ClinGen dosage sensitivity table: {r}")
    lines = r.content.decode('UTF-8').split("\n")
    header_idx = None
    for i, line in enumerate(lines):
        if '"GENE SYMBOL"' in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find header line containing 'GENE SYMBOL' in ClinGen dosage download")
    data_start = header_idx + 1
    if data_start < len(lines) and lines[data_start].startswith('"+++'):
        data_start += 1
    table_contents = "\n".join([lines[header_idx]] + lines[data_start:])
    df = pd.read_csv(StringIO(table_contents))
    df = df[["HGNC ID", "HAPLOINSUFFICIENCY"]]
    df = df[~df["HAPLOINSUFFICIENCY"].isin([
        "No Evidence for Haploinsufficiency", 
        "Dosage Sensitivity Unlikely",
        "Little Evidence for Haploinsufficiency",
    ])]
    return df


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_clingen_haploinsufficient_genes_table()
    print("Haploinsufficient Genes Table columns:")
    print(df.iloc[0])

    unknown_hgnc_ids = set(df["HGNC ID"]) - set(get_hgnc_to_ensg_id_map().keys())
    assert len(unknown_hgnc_ids) == 0, f"Unknown HGNC ids in haploinsufficient genes table: {', '.join(unknown_hgnc_ids)}"


    df = get_clingen_gene_disease_validity_table()

    unknown_hgnc_ids = set(df["GENE ID (HGNC)"]) - set(get_hgnc_to_ensg_id_map().keys())
    assert len(unknown_hgnc_ids) == 0, f"Unknown HGNC ids in gene-disease validity table: {', '.join(unknown_hgnc_ids)}"


    print("Gene-Disease Validity Table columns:")
    print(df.iloc[0])

    print(f"Inheritance column:")
    print(df["MOI"].value_counts())




"""
5/19/2025 ClinGen Gene-Disease Validity Table:

Columns:
0: GENE SYMBOL
1: GENE ID (HGNC)
2: DISEASE LABEL
3: DISEASE ID (MONDO)
4: MOI
5: SOP
6: CLASSIFICATION
7: ONLINE REPORT
8: CLASSIFICATION DATE
9: GCEP

Example:

GENE SYMBOL                                                        AARS1
GENE ID (HGNC)                                                   HGNC:20
DISEASE LABEL                 Charcot-Marie-Tooth disease axonal type 2N
DISEASE ID (MONDO)                                         MONDO:0013212
MOI                                                                   AD
SOP                                                                SOP10
CLASSIFICATION                                                Definitive
ONLINE REPORT          https://search.clinicalgenome.org/kb/gene-vali...
CLASSIFICATION DATE                             2024-03-14T16:00:00.000Z
GCEP                   Charcot-Marie-Tooth Disease Gene Curation Expe...


Unique GeneIds: 2,564


CLASSIFICATION
Definitive                       1932
Limited                           461
Moderate                          367
Disputed                          177
Strong                             70
No Known Disease Relationship      48
Refuted                            40


MOI
AR    1519
AD    1260
XL     182
MT      54
SD      51
UD      29
"""
