import pandas as pd
import requests
from io import StringIO

from annotation_utils.cache_utils import cache_data_table

URL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=gd_mgd_id&col=gd_pubmed_ids&col=gd_locus_type&col=gd_pub_eg_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

@cache_data_table
def _download_hgnc_table(url):
    r = requests.get(url)
    if not r.ok:
        raise Exception(f"Failed to download {url}: {r}")

    table_contents = StringIO(r.content.decode('UTF-8'))
    return pd.read_table(table_contents)


def get_hgnc_table():
    """Download the HGNC table from https://www.genenames.org/download/custom/ and return it as a pandas DataFrame.

    The URL is passed as an argument to the cached downloader so it becomes part of the
    @cache_data_table cache key. Changing the requested columns (e.g. adding NCBI Gene ID) then
    invalidates any stale cache instead of silently reusing a table that lacks the new column.
    """
    return _download_hgnc_table(URL)

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
        # When an Ensembl gene ID maps to multiple HGNC IDs, keep only the first one. The
        # downstream consumer interpolates this value into single-ID URLs (genenames.org,
        # ClinGen, etc.) on the per-gene page, which break for comma-joined multi-IDs.
        df_hgnc = df_hgnc.groupby("Ensembl gene ID", as_index=False).agg({"HGNC ID": "first"})

    ENSG_to_HGNC_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["HGNC ID"]))

    return ENSG_to_HGNC_map


def _map_each_listed_symbol_to_ensg_ids(df_hgnc, symbols_column):
    """Return {symbol: set of ENSG ids} for the comma-separated symbols in one column of the HGNC table."""
    symbol_to_ensg_ids = {}
    for ensg_id, symbols in zip(df_hgnc["Ensembl gene ID"], df_hgnc[symbols_column]):
        if not isinstance(symbols, str):
            continue
        for symbol in symbols.split(","):
            if symbol.strip():
                symbol_to_ensg_ids.setdefault(symbol.strip(), set()).add(ensg_id)
    return symbol_to_ensg_ids


def build_gene_symbol_to_ensg_id_map(df_hgnc):
    """Map gene symbols, including outdated ones, to Ensembl gene ids.

    For sources that identify genes only by symbol, some of which may be out of date. A symbol resolves,
    in order of preference, to the gene that currently has it as its approved symbol, else to the one gene
    that lists it as a previous symbol, else to the one gene that lists it as an alias. Symbols that are
    previous symbols or aliases of more than one gene are left out because they are ambiguous.

    Args:
        df_hgnc (pandas.DataFrame): the HGNC table returned by get_hgnc_table()

    Returns:
        dict: {gene symbol: Ensembl gene id}
    """
    df_hgnc = df_hgnc[df_hgnc["Ensembl gene ID"].notna()]
    previous_symbol_to_ensg_ids = _map_each_listed_symbol_to_ensg_ids(df_hgnc, "Previous symbols")
    alias_to_ensg_ids = _map_each_listed_symbol_to_ensg_ids(df_hgnc, "Alias symbols")

    symbol_to_ensg_id = {
        alias: ensg_ids.pop() for alias, ensg_ids in alias_to_ensg_ids.items()
        if len(ensg_ids) == 1 and alias not in previous_symbol_to_ensg_ids
    }
    symbol_to_ensg_id.update({
        previous_symbol: ensg_ids.pop() for previous_symbol, ensg_ids in previous_symbol_to_ensg_ids.items()
        if len(ensg_ids) == 1
    })
    # Approved symbols are applied last so that they override any gene that lists the same symbol as a
    # previous symbol or alias. When an approved symbol appears on more than one row, the first row wins.
    for approved_symbol, ensg_id in zip(df_hgnc["Approved symbol"][::-1], df_hgnc["Ensembl gene ID"][::-1]):
        symbol_to_ensg_id[approved_symbol] = ensg_id

    return symbol_to_ensg_id


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
13: NCBI Gene ID

(44109, 14)

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