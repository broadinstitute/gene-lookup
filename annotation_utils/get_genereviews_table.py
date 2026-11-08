"""Download the list of GeneReviews chapters for each gene from NCBI.

GeneReviews (https://www.ncbi.nlm.nih.gov/books/NBK1116/) are expert-written clinical reviews of genetic
conditions, each published as an NCBI Bookshelf chapter with an NBK id (eg. NBK1103 is
https://www.ncbi.nlm.nih.gov/books/NBK1103/).

Example rows of the two files used here (they are Latin-1 encoded, eg. "Strømme Syndrome"):

    NBKid_shortname_genesymbol.txt:   #NBK_id  GR_shortname       genesymbol
                                      NBK1103  trimethylaminuria  FMO3
    GRtitle_shortname_NBKid.txt:      #GR_shortname      GR_Title                   NBK_id   PMID
                                      trimethylaminuria  Primary Trimethylaminuria  NBK1103  20301282
"""

import io

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise, utc_today

GENEREVIEWS_GENE_SYMBOLS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/NBKid_shortname_genesymbol.txt"
GENEREVIEWS_TITLES_URL = "https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt"


def combine_genereviews_gene_symbols_and_titles(df_gene_symbols, df_titles):
    """Join each gene's GeneReviews chapters to their titles.

    Args:
        df_gene_symbols (pandas.DataFrame): NBKid_shortname_genesymbol.txt
        df_titles (pandas.DataFrame): GRtitle_shortname_NBKid.txt

    Returns:
        pandas.DataFrame: one row per gene and chapter, with columns gene_symbol, nbk_id, title. Rows whose
            chapter has no NBK id, or whose gene symbol is "Not applicable" (chapters about conditions not
            tied to specific genes), are dropped.
    """
    df = df_gene_symbols.rename(columns={"#NBK_id": "nbk_id", "genesymbol": "gene_symbol"})
    df = df[df["nbk_id"].notna() & df["gene_symbol"].notna() & (df["gene_symbol"] != "Not applicable")]

    nbk_id_to_title = dict(zip(df_titles["NBK_id"], df_titles["GR_Title"]))
    missing_titles = set(df["nbk_id"]) - set(nbk_id_to_title)
    if missing_titles:
        raise ValueError(f"GeneReviews chapters without a title: {sorted(missing_titles)}")

    df = df.assign(title=df["nbk_id"].map(nbk_id_to_title))
    return df[["gene_symbol", "nbk_id", "title"]].drop_duplicates().reset_index(drop=True)


def _download_latin1_table(url):
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return pd.read_table(io.BytesIO(response.content), encoding="latin-1")


@cache_data_table
def get_genereviews_table():
    """Download the GeneReviews chapters for each gene.

    Returns:
        pandas.DataFrame: see combine_genereviews_gene_symbols_and_titles. A gene can have several chapters.
    """
    try:
        df = combine_genereviews_gene_symbols_and_titles(
            _download_latin1_table(GENEREVIEWS_GENE_SYMBOLS_URL), _download_latin1_table(GENEREVIEWS_TITLES_URL))
    except Exception as e:
        return return_stale_cached_table_or_raise("get_genereviews_table", "GeneReviews", e)

    # These files carry no release date, so the download date is all we have.
    record_source_version("genereviews", date_label=f"downloaded: {utc_today()}")

    return df


if __name__ == "__main__":
    print(get_genereviews_table())
