import pandas as pd
import requests
from io import StringIO

from annotation_utils.cache_utils import cache_data_table

URL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=gd_mgd_id&col=gd_pubmed_ids&col=gd_locus_type&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

@cache_data_table
def get_hgnc_table():
    """Download the HGNC table from https://www.genenames.org/download/custom/ and return it as a pandas DataFrame"""
    r = requests.get(URL)
    if not r.ok:
        raise Exception(f"Failed to download {URL}: {r}")

    table_contents = StringIO(r.content.decode('UTF-8'))
    return pd.read_table(table_contents)

def _remove_duplicate_ensg_ids(df_hgnc):
    before = len(df_hgnc)
    locus_types_of_genes_with_duplicate_ENSG_ids = df_hgnc[df_hgnc["Ensembl gene ID"].duplicated()]["Locus type"].value_counts()
    df_hgnc = df_hgnc[~df_hgnc["Ensembl gene ID"].duplicated()]
    print(f"Discarded {before - len(df_hgnc):,d} out of {before:,d} ({(before - len(df_hgnc))/before*100:.1f}%) which had duplicate Ensembl gene IDs. Their locus types were:")
    print(locus_types_of_genes_with_duplicate_ENSG_ids)
    return df_hgnc

def get_hgnc_to_ensg_id_map():
    # there is exactly one Ensembl gene ID per HGNC ID, so we can just use the HGNC ID as the key
    df_hgnc = get_hgnc_table()
    n_duplicates = df_hgnc["HGNC ID"].duplicated().sum()
    assert n_duplicates == 0, f"Found {n_duplicates:,d} duplicate HGNC IDs in the HGNC table: " + ", ".join(df_hgnc[df_hgnc["HGNC ID"].duplicated()]["HGNC ID"].unique())
    HGNC_to_ENSG_map = dict(zip(df_hgnc["HGNC ID"], df_hgnc["Ensembl gene ID"]))

    return HGNC_to_ENSG_map


def get_ensg_id_to_hgnc_id_map(remove_duplicate_ensg_ids=False):
    df_hgnc = get_hgnc_table()
    if remove_duplicate_ensg_ids:
        # this is a one to many mapping, so we need to remove the duplicates
        df_hgnc = _remove_duplicate_ensg_ids(df_hgnc)
    else:
        # group by Ensembl gene ID and join the HGNC IDs by comma
        df_hgnc = df_hgnc.groupby("Ensembl gene ID").agg({"HGNC ID": lambda x: ",".join(x)}).reset_index()

    ENSG_to_HGNC_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["HGNC ID"]))

    return ENSG_to_HGNC_map

if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_hgnc_table()
    print(df)


"""
5/19/2025HGNC Table:

0: HGNC ID
1: Approved symbol
2: Approved name
3: Status
4: Previous symbols
5: Alias symbols
6: Chromosome
7: Accession numbers
8: RefSeq IDs
9: Ensembl gene ID
10: Mouse genome database ID
11: Pubmed IDs
12: Locus type

(44109, 13)

Example:

HGNC ID                                        HGNC:5
Approved symbol                                  A1BG
Approved name                  alpha-1-B glycoprotein
Status                                       Approved
Previous symbols                                  NaN
Alias symbols                                     NaN
Chromosome                                   19q13.43
Accession numbers                                 NaN
RefSeq IDs                                  NM_130786
Ensembl gene ID                       ENSG00000121410
Mouse genome database ID                  MGI:2152878
Pubmed IDs                                    2591067
Locus type                  gene with protein product

"""