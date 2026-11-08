"""Download ClinGen's clinical actionability curations for adult and pediatric contexts.

ClinGen's Actionability Working Group (https://actionability.clinicalgenome.org/ac/) scores how actionable a
gene-condition pair is: for each outcome-intervention pair it scores severity, likelihood, effectiveness and
nature of the intervention, and sums them into an overall score from 0 to 12, followed by letters grading the
evidence for likelihood and effectiveness (eg. "10CB"). Curations are done separately for adults and children.

Example row of the flat summary for the Adult context (fields not used here are omitted):

    docId           AC1084
    context         Adult
    curationType    Gene-Condition
    geneOrVariant   CYP27A1
    disease         Cerebrotendinous xanthomatosis
    status-overall  Released
    overall         10CN

A report can cover several genes (geneOrVariant "MLH1,MSH2,MSH6,PMS2") and has one row per outcome-intervention
pair. Each report is at https://actionability.clinicalgenome.org/ac/<context>/ui/stg2SummaryRpt?doc=<docId>
"""

import re

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise, utc_today

ACTIONABILITY_CONTEXTS = ["Adult", "Pediatric"]
ACTIONABILITY_SUMMARY_URL_TEMPLATE = "https://actionability.clinicalgenome.org/ac/{context}/api/summ?flavor=flat"


def parse_overall_actionability_score(overall_score_text):
    """Return the numeric part of an overall actionability score (10 for "10CB"), or None if it has none."""
    match = re.match(r"^(\d+)", str(overall_score_text))
    return int(match.group(1)) if match else None


def convert_actionability_summaries_to_gene_table(df_summary_rows):
    """Turn the flat actionability summary rows into one row per released gene-condition report and gene.

    Args:
        df_summary_rows (pandas.DataFrame): rows of the flat summaries for one or more contexts, one row per
            report and outcome-intervention pair

    Returns:
        pandas.DataFrame: columns gene_symbol, context, report_id, disease, max_overall_score (the highest
            overall score among the report's outcome-intervention pairs, or missing if none is scored)
    """
    df = df_summary_rows[
        (df_summary_rows["status-overall"] == "Released") & (df_summary_rows["curationType"] == "Gene-Condition")
    ]
    df = df.assign(
        overall_score=df["overall"].map(parse_overall_actionability_score),
        gene_symbol=df["geneOrVariant"].str.split(","),
    ).explode("gene_symbol")
    df = df.assign(gene_symbol=df["gene_symbol"].str.strip())
    # A report without genes lists the placeholder "(ERR: No gene(s) listed for topic?)" instead of a symbol
    df = df[(df["gene_symbol"] != "") & ~df["gene_symbol"].str.startswith("(")]

    return df.groupby(["gene_symbol", "context", "docId"], as_index=False).agg(
        disease=("disease", "first"),
        max_overall_score=("overall_score", "max"),
    ).rename(columns={"docId": "report_id"})[["gene_symbol", "context", "report_id", "disease", "max_overall_score"]]


@cache_data_table
def get_clingen_actionability_table():
    """Download the released ClinGen actionability curations for adults and children.

    Returns:
        pandas.DataFrame: see convert_actionability_summaries_to_gene_table
    """
    try:
        summary_tables = []
        for context in ACTIONABILITY_CONTEXTS:
            response = requests.get(ACTIONABILITY_SUMMARY_URL_TEMPLATE.format(context=context), timeout=120)
            response.raise_for_status()
            summary = response.json()
            summary_tables.append(pd.DataFrame(summary["rows"], columns=summary["columns"]))
        df = convert_actionability_summaries_to_gene_table(pd.concat(summary_tables, ignore_index=True))
    except Exception as e:
        return return_stale_cached_table_or_raise("get_clingen_actionability_table", "ClinGen actionability", e)

    # The summaries carry per-report release dates but no overall release, so the download date is used.
    record_source_version("clingen_actionability", date_label=f"downloaded: {utc_today()}")

    return df


if __name__ == "__main__":
    print(get_clingen_actionability_table())
