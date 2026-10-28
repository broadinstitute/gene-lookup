from annotation_utils.cache_utils import cache_data_table, read_cached_table, record_source_version, utc_today
import re
import requests
import pandas as pd

"""
Example record from https://rpldb.org/intolerome/term/geneData (the "_source" of each element of hits.hits):

{'gene': {'id': '6299', 'name': 'SALL1'},
 'omim_gene': '602218',
 'gene_effect': ['Fetal'],
 'variants': ['CNV', 'SNV'],
 'mode_of_function': 'LoF',
 'coordinates': {'id': 'chr16', 'data': '51,169,893-51,185,181'},
 'disease': 'Townes-Brocks branchiootorenal-like syndrome, 107480 (3), Autosomal dominant, Townes-Brocks syndrome 1, 107480 (3), Autosomal dominant',
 'omim_phenotype': ['107480'],
 'inheritance': ['AD'],
 'prenatal_system_affected': ['MULTIPLE >3'],
 'prenatal_phenotypes': ['Neonatal demise.'],
 'fetal_sex': ['M', 'F'],
 'gene_category_and_human_evidence_rating': 'Known',
 'pmid': ['30993004', '11484202'],
 '1_trimester': '',
 '2_trimester': '',
 '3_trimester': '',
 'postnatal_lethality': 'x'}
"""

INTOLEROME_GENE_DATA_URL = "https://rpldb.org/intolerome/term/geneData"

# Abbreviate Intolerome's X-linked codes the same way get_omim_table.py and get_gencc_table.py do, so the
# combined inheritance column doesn't list the same mode twice (eg. Intolerome's "XLR" next to OMIM's "XR")
INHERITANCE_ABBREVIATIONS = {
    "XLR": "XR",
    "XLD": "XD",
    "X-linked": "XR",
}

# Intolerome marks each period in which deaths were reported with an "x" in that period's field
LETHALITY_TIMING_FIELDS_AND_LABELS = [
    ("1_trimester", "1st trimester"),
    ("2_trimester", "2nd trimester"),
    ("3_trimester", "3rd trimester"),
    ("postnatal_lethality", "postnatal"),
]


def _join_non_empty_values(values, separator=", "):
    return separator.join(v.strip() for v in values if v and v.strip())


def _clean_free_text(text):
    """Replace the literal "\\n" sequences and runs of whitespace in Intolerome's free text with single spaces."""
    return re.sub(r"\s+", " ", text.replace("\\n", " ")).strip()


def _convert_intolerome_record_to_row(record):
    return {
        "ncbi_gene_id": record["gene"]["id"],
        "gene_symbol": record["gene"]["name"],
        "evidence_category": record["gene_category_and_human_evidence_rating"],
        "lethality_timing": _join_non_empty_values(
            [label for field, label in LETHALITY_TIMING_FIELDS_AND_LABELS if record[field].strip()]),
        "gene_effect": _join_non_empty_values(record["gene_effect"]),
        # Joined with "; " so that summarize_inheritance in generate_combined_gene_table.py splits a record
        # listing several modes (eg. AD and AR) into separate codes
        "inheritance": _join_non_empty_values(
            [INHERITANCE_ABBREVIATIONS.get(i.strip(), i) for i in record["inheritance"]], separator="; "),
        "disease": _clean_free_text(record["disease"]),
        "omim_phenotype_id": _join_non_empty_values(record["omim_phenotype"]),
        "prenatal_phenotypes": _clean_free_text(_join_non_empty_values(record["prenatal_phenotypes"], separator=" ")),
        "system_affected": _join_non_empty_values(record["prenatal_system_affected"]),
        "pmids": _join_non_empty_values(record["pmid"]),
    }


@cache_data_table
def get_intolerome_table():
    """Download the Intolerome list of genes linked to pregnancy loss, stillbirth, and neonatal death.

    Source: Yatsenko et al. 2026 (https://pubmed.ncbi.nlm.nih.gov/41793222/), database at
    https://rpldb.org/intolerome/

    Returns:
        pandas.DataFrame: one row per Intolerome record, identified by NCBI gene id. A few genes have
            two records (eg. one "Known" and one "Emerging"), so ncbi_gene_id is not unique.
    """
    try:
        response = requests.get(INTOLEROME_GENE_DATA_URL, headers={"accept": "application/json"}, timeout=120)
        if not response.ok:
            raise ValueError(f"Failed to get data from {INTOLEROME_GENE_DATA_URL}: HTTP {response.status_code} {response.text[:200]}")

        hits = response.json()["hits"]
        records = [hit["_source"] for hit in hits["hits"]]
        if not records:
            raise ValueError(f"Failed to get records from {INTOLEROME_GENE_DATA_URL}: {hits}")

        # The endpoint runs an Elasticsearch query, which returns at most one page of results. Fail rather
        # than load a partial gene list if the database ever outgrows that page.
        if len(records) != hits["total"]["value"]:
            raise ValueError(f"{INTOLEROME_GENE_DATA_URL} returned only {len(records):,d} out of {hits['total']['value']:,d} records")

        output_records = [_convert_intolerome_record_to_row(record) for record in records]
    except Exception as e:
        # This is an internal endpoint of the Intolerome website rather than a supported download, so its
        # format may change without notice. Fall back to the last successfully cached Intolerome table
        # (even if it's stale) rather than failing the whole pipeline.
        df_stale = read_cached_table("get_intolerome_table")
        if df_stale is not None:
            print(f"WARNING: Intolerome download failed ({e!r}); reusing previously cached Intolerome table ({len(df_stale):,d} rows)")
            return df_stale
        raise

    # This endpoint returns no version or build date, so the download date is all we have.
    record_source_version("intolerome", date_label=f"downloaded: {utc_today()}")

    return pd.DataFrame(output_records)


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_intolerome_table()
    print(df)
