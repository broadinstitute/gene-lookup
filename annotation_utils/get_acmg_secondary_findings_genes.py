"""Download the ACMG secondary findings gene list from ClinGen.

The American College of Medical Genetics and Genomics (ACMG) secondary findings (SF) list names the genes
whose pathogenic variants should be reported from clinical exome and genome sequencing even when they are
unrelated to the reason for testing. ClinGen hosts the current list at
https://search.clinicalgenome.org/kb/genes/acmgsf, which loads its rows from the JSON endpoint used below.

Example row from the endpoint (fields not used here are omitted):

    {"symbol": "ABCD1", "hgnc_id": "HGNC:61", ...}
"""

import re
import time

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise

CLINGEN_ACMG_SECONDARY_FINDINGS_PAGE_URL = "https://search.clinicalgenome.org/kb/genes/acmgsf"
CLINGEN_ACMG_SECONDARY_FINDINGS_GENES_URL = "https://search.clinicalgenome.org/api/genes/acmg?offset=0&limit=1000&order=asc"

# In November 2026 the endpoint answered most requests with HTTP 500 and roughly one in four with the data,
# regardless of the request headers, so it is retried before falling back to the cached table.
MAX_DOWNLOAD_ATTEMPTS = 20
SECONDS_BETWEEN_DOWNLOAD_ATTEMPTS = 5


def _get_with_retries(url, headers=None):
    for attempt in range(1, MAX_DOWNLOAD_ATTEMPTS + 1):
        try:
            response = requests.get(url, headers=headers, timeout=60)
            response.raise_for_status()
            return response
        except Exception as e:
            if attempt == MAX_DOWNLOAD_ATTEMPTS:
                raise
            print(f"  Attempt {attempt} to download {url} failed ({e!r}). Retrying in {SECONDS_BETWEEN_DOWNLOAD_ATTEMPTS} seconds...")
            time.sleep(SECONDS_BETWEEN_DOWNLOAD_ATTEMPTS)


def parse_acmg_secondary_findings_list_version(page_html):
    """Return the ACMG SF list version (eg. "3.3") named in the ClinGen page's "ACMG SF v3.3 Curations"
    heading, or "" if the page doesn't name one."""
    versions = re.findall(r"ACMG SF v(\d+(?:\.\d+)+) Curations", page_html)
    return max(versions, key=lambda v: tuple(int(part) for part in v.split("."))) if versions else ""


def convert_acmg_secondary_findings_response_to_table(response_json):
    """Convert the ClinGen endpoint's JSON response to a table of genes.

    Args:
        response_json (dict): the endpoint's response, with "total" and "rows" keys

    Returns:
        pandas.DataFrame: one row per gene, with columns hgnc_gene_id and gene_symbol
    """
    rows = response_json["rows"]
    if len(rows) != response_json["total"]:
        raise ValueError(f"The ClinGen ACMG SF endpoint returned {len(rows)} out of {response_json['total']} genes")

    return pd.DataFrame({
        "hgnc_gene_id": [row["hgnc_id"] for row in rows],
        "gene_symbol": [row["symbol"] for row in rows],
    }).drop_duplicates().reset_index(drop=True)


@cache_data_table
def get_acmg_secondary_findings_genes():
    """Download the genes on the current ACMG secondary findings list.

    Returns:
        pandas.DataFrame: one row per gene, with columns hgnc_gene_id and gene_symbol
    """
    try:
        response = _get_with_retries(CLINGEN_ACMG_SECONDARY_FINDINGS_GENES_URL, headers={"Accept": "application/json"})
        df = convert_acmg_secondary_findings_response_to_table(response.json())
    except Exception as e:
        return return_stale_cached_table_or_raise("get_acmg_secondary_findings_genes", "ACMG secondary findings", e)

    # The version only labels the website heading, so a failure to read it leaves the previous label in place
    # rather than failing the download.
    try:
        list_version = parse_acmg_secondary_findings_list_version(_get_with_retries(CLINGEN_ACMG_SECONDARY_FINDINGS_PAGE_URL).text)
    except Exception as e:
        print(f"WARNING: could not read the ACMG secondary findings list version from {CLINGEN_ACMG_SECONDARY_FINDINGS_PAGE_URL} ({e!r})")
        list_version = ""
    if list_version:
        record_source_version("acmg_secondary_findings", version=f"v{list_version}")

    return df


if __name__ == "__main__":
    print(get_acmg_secondary_findings_genes())
