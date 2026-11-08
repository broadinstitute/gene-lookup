"""Download the gene burden evidence from the latest Open Targets Platform release.

Gene burden evidence (https://platform-docs.opentargets.org/evidence#gene-burden) is a set of rare-variant
burden test associations between genes and traits, mostly from UK Biobank exome sequencing analyses
(AstraZeneca PheWAS Portal, Genebass, Regeneron, BRaVa), plus disease-specific studies such as the Autism
Sequencing Consortium and SCHEMA. Open Targets maps each trait to an EFO, MONDO, HP or other ontology id.

Example record (columns not used here are omitted):

    targetId                   ENSG00000118137
    diseaseId                  EFO_0022282
    projectId                  AstraZeneca PheWAS Portal
    cohortId                   UK Biobank 470k
    literature                 ['34375979']
    pValueMantissa             1.677
    pValueExponent             -15
    directionOnTrait           protect

About 200 records from individual publications have no projectId and name only their cohort (eg. "INTERVAL",
"SPARK, ASC, MSSNG, SSC") or only their PubMed ID.
"""

import io
import re

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise

OPEN_TARGETS_PLATFORM_RELEASES_URL = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/"

# Passed to get_open_targets_gene_burden_table so that it becomes part of the cache key. Bump it whenever the
# columns returned by convert_gene_burden_evidence_to_table change, so that a table cached with the old columns
# is downloaded again instead of reused. Version 2 added study_pubmed_id.
OPEN_TARGETS_GENE_BURDEN_TABLE_COLUMNS_VERSION = 2


def find_latest_open_targets_release(directory_listing_html):
    """Return the newest release (eg. "26.09") in the Open Targets Platform FTP directory listing.

    Args:
        directory_listing_html (str): HTML of the OPEN_TARGETS_PLATFORM_RELEASES_URL directory listing

    Returns:
        str: the release name in YY.MM form
    """
    releases = re.findall(r'href="(\d{2}\.\d{2})/"', directory_listing_html)
    if not releases:
        raise ValueError(f"No releases found in the {OPEN_TARGETS_PLATFORM_RELEASES_URL} listing")
    return max(releases, key=lambda release: tuple(int(part) for part in release.split(".")))


def find_parquet_file_names(directory_listing_html):
    """Return the sorted names of the .parquet files in an FTP directory listing."""
    return sorted(set(re.findall(r'href="([^"/]+\.parquet)"', directory_listing_html)))


def _get_first_pubmed_id(pubmed_ids):
    """Return the first PubMed ID in the evidence record's literature list, or "" if the list is empty."""
    if pubmed_ids is not None and len(pubmed_ids) > 0:
        return str(pubmed_ids[0])
    return ""


def _name_study(project_id, cohort_id, pubmed_ids):
    """Return the project name, else the cohort name, else "PMID:<first PubMed ID>", else ""."""
    for name in (project_id, cohort_id):
        if isinstance(name, str) and name.strip():
            return " ".join(name.split())
    if _get_first_pubmed_id(pubmed_ids):
        return f"PMID:{_get_first_pubmed_id(pubmed_ids)}"
    return ""


def convert_gene_burden_evidence_to_table(df_evidence):
    """Select and normalize the gene burden evidence columns used by gene-lookup.

    Args:
        df_evidence (pandas.DataFrame): the evidence_gene_burden parquet table

    Returns:
        pandas.DataFrame: one row per association, with columns gene_id, disease_id (with ":" in place of
            Open Targets' "_", eg. "MONDO:0005258"), study, study_pubmed_id (the PubMed ID of the study's
            paper, or "" if Open Targets lists none), p_value, direction_on_trait ("risk", "protect", or
            missing)
    """
    return pd.DataFrame({
        "gene_id": df_evidence["targetId"],
        "disease_id": df_evidence["diseaseId"].str.replace("_", ":", n=1, regex=False),
        "study": [_name_study(project_id, cohort_id, pubmed_ids) for project_id, cohort_id, pubmed_ids in zip(
            df_evidence["projectId"], df_evidence["cohortId"], df_evidence["literature"])],
        "study_pubmed_id": [_get_first_pubmed_id(pubmed_ids) for pubmed_ids in df_evidence["literature"]],
        # p-values below ~1e-308 underflow to 0.0, which still sorts as the most significant
        "p_value": df_evidence["pValueMantissa"] * 10.0 ** df_evidence["pValueExponent"],
        "direction_on_trait": df_evidence["directionOnTrait"],
    })


def _get_url_text(url):
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


@cache_data_table
def get_open_targets_gene_burden_table(table_columns_version):
    """Download the gene burden evidence from the latest Open Targets Platform release.

    Args:
        table_columns_version (int): OPEN_TARGETS_GENE_BURDEN_TABLE_COLUMNS_VERSION. It is not used here, but
            cache_data_table includes it in the cache key, so bumping it makes the pipeline download a fresh table.

    Returns:
        pandas.DataFrame: one row per gene-trait association (see convert_gene_burden_evidence_to_table)
    """
    try:
        release = find_latest_open_targets_release(_get_url_text(OPEN_TARGETS_PLATFORM_RELEASES_URL))
        evidence_directory_url = f"{OPEN_TARGETS_PLATFORM_RELEASES_URL}{release}/output/evidence_gene_burden/"
        parquet_file_names = find_parquet_file_names(_get_url_text(evidence_directory_url))
        if not parquet_file_names:
            raise ValueError(f"No parquet files found in {evidence_directory_url}")

        evidence_tables = []
        for parquet_file_name in parquet_file_names:
            print(f"Downloading {evidence_directory_url}{parquet_file_name}")
            response = requests.get(f"{evidence_directory_url}{parquet_file_name}", timeout=300)
            response.raise_for_status()
            evidence_tables.append(pd.read_parquet(io.BytesIO(response.content)))
        df = convert_gene_burden_evidence_to_table(pd.concat(evidence_tables, ignore_index=True))
    except Exception as e:
        return return_stale_cached_table_or_raise("get_open_targets_gene_burden_table", "Open Targets gene burden", e)

    record_source_version("open_targets_gene_burden", date_label=f"release: {release}")

    return df


if __name__ == "__main__":
    print(get_open_targets_gene_burden_table(OPEN_TARGETS_GENE_BURDEN_TABLE_COLUMNS_VERSION))
