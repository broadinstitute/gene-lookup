"""Download Gene2Phenotype (G2P) gene-disease records from EBI's monthly data releases.

G2P (https://www.ebi.ac.uk/gene2phenotype/) curates gene-disease associations across its DD, Eye, Skin,
Skeletal, Cardiac, Cancer and Ear panels. Beyond the gene-disease pairs, each record states the molecular
mechanism (loss of function, gain of function, dominant negative) and the variant consequence.

Example record from allG2P_2026-10-28.csv.gz (columns not used here are omitted):

    g2p id                  G2P00001
    gene symbol             HMX1
    hgnc id                 5017
    disease name            HMX1-related oculoauricular syndrome
    disease mim             612109.0
    disease MONDO           MONDO:0012802
    allelic requirement     biallelic_autosomal
    confidence              definitive
    variant consequence     absent gene product
    molecular mechanism     loss of function
    panel                   DD; Eye
"""

import io
import re

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise

G2P_DATA_DOWNLOADS_URL = "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/"

# Map G2P allelic requirements to the inheritance codes the other sources use (see get_gencc_table.py), so
# the combined inheritance column doesn't list the same mode twice under different names.
G2P_ALLELIC_REQUIREMENT_TO_INHERITANCE = {
    "biallelic_autosomal": "AR",
    "monoallelic_autosomal": "AD",
    "monoallelic_X_hemizygous": "XR",
    "monoallelic_X_heterozygous": "XD",
    "monoallelic_X": "XL",
    "monoallelic_Y_hemizygous": "YL",
    "mitochondrial": "MITO",
    "biallelic_PAR": "Pseudoautosomal recessive",
    "monoallelic_PAR": "Pseudoautosomal dominant",
}


# G2P abbreviates its developmental disorders panel as "DD"; gene-lookup stores and shows the full name instead.
G2P_PANEL_ABBREVIATION_TO_NAME = {"DD": "Developmental disorders"}


def spell_out_abbreviated_g2p_panel_names(panel_values):
    """Replace each abbreviated panel name in G2P_PANEL_ABBREVIATION_TO_NAME with its full name.

    This runs on the table returned by get_g2p_table, so it also applies to a table read from the cache.

    Args:
        panel_values (pandas.Series): the panel column, with a record's panels joined with ", ", eg. "DD, Eye"

    Returns:
        pandas.Series: the same values with abbreviations spelled out, eg. "Developmental disorders, Eye"
    """
    # ANALYSIS_OK[imputation]: "" marks a record with no panel, as normalize_nulls does when records are aggregated
    return panel_values.fillna("").astype(str).map(
        lambda value: ", ".join(G2P_PANEL_ABBREVIATION_TO_NAME.get(panel, panel) for panel in value.split(", ")))


def find_latest_g2p_release_date(directory_listing_html):
    """Return the newest release directory name (eg. "2026_10_28") in the G2P data downloads listing.

    Args:
        directory_listing_html (str): HTML of the G2P_DATA_DOWNLOADS_URL directory listing

    Returns:
        str: the release date in the YYYY_MM_DD form used for the directory names
    """
    release_dates = re.findall(r'href="(\d{4}_\d{2}_\d{2})/"', directory_listing_html)
    if not release_dates:
        raise ValueError(f"No release directories found in the {G2P_DATA_DOWNLOADS_URL} listing")
    return max(release_dates)


def convert_g2p_records_to_gene_disease_table(df_g2p):
    """Select and normalize the G2P columns used by gene-lookup.

    Args:
        df_g2p (pandas.DataFrame): the allG2P_<date>.csv table, one row per G2P record

    Returns:
        pandas.DataFrame: one row per G2P record, with HGNC ids in "HGNC:<number>" form, inheritance codes,
            and multi-valued fields joined with ", " so that the "; " separator used when aggregating
            records per gene keeps one token per record.
    """
    unmapped_allelic_requirements = set(df_g2p["allelic requirement"].dropna()) - set(G2P_ALLELIC_REQUIREMENT_TO_INHERITANCE)
    if unmapped_allelic_requirements:
        print(f"WARNING: G2P allelic requirements without an inheritance code will be passed through unchanged: {sorted(unmapped_allelic_requirements)}")

    def join_multiple_values_with_commas(values):
        # ANALYSIS_OK[imputation]: "" marks a missing value; without it astype(str) would write the text "nan"
        return values.fillna("").astype(str).str.replace(r"\s*;\s*", ", ", regex=True)

    return pd.DataFrame({
        "hgnc_gene_id": "HGNC:" + df_g2p["hgnc id"].astype(int).astype(str),
        "g2p_id": df_g2p["g2p id"],
        "disease_name": join_multiple_values_with_commas(df_g2p["disease name"]),
        "inheritance": df_g2p["allelic requirement"].replace(G2P_ALLELIC_REQUIREMENT_TO_INHERITANCE),
        "confidence": df_g2p["confidence"],
        "molecular_mechanism": join_multiple_values_with_commas(df_g2p["molecular mechanism"]),
        "variant_consequence": join_multiple_values_with_commas(df_g2p["variant consequence"]),
        "panel": join_multiple_values_with_commas(df_g2p["panel"]),
    })


@cache_data_table
def get_g2p_table():
    """Download the latest monthly G2P release covering all panels.

    Returns:
        pandas.DataFrame: one row per G2P record (see convert_g2p_records_to_gene_disease_table). A gene can
            have several records, one per disease or allelic requirement.
    """
    try:
        listing_response = requests.get(G2P_DATA_DOWNLOADS_URL, timeout=60)
        listing_response.raise_for_status()
        release_date = find_latest_g2p_release_date(listing_response.text)
        url = f"{G2P_DATA_DOWNLOADS_URL}{release_date}/allG2P_{release_date.replace('_', '-')}.csv.gz"
        print(f"Downloading {url}")
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        df_g2p = pd.read_csv(io.BytesIO(response.content), compression="gzip")
    except Exception as e:
        return return_stale_cached_table_or_raise("get_g2p_table", "G2P", e)

    record_source_version("g2p", date_label=f"released: {release_date.replace('_', '-')}")

    return convert_g2p_records_to_gene_disease_table(df_g2p)


if __name__ == "__main__":
    print(get_g2p_table())
