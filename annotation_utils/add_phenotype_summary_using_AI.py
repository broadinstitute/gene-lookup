import argparse
import pandas as pd
from dotenv import load_dotenv
import sys
load_dotenv()

from llm_utils.text_completion import ask_gemini

parser = argparse.ArgumentParser()
parser.add_argument("--output-path", type=str, help="The output table path")
parser.add_argument("combined_table", type=str, help="The combined TSV table of genes and their phenotypes")
args = parser.parse_args()

if not args.output_path:
    args.output_path = args.combined_table.replace(".tsv", "_with_phenotype_summary.tsv")

prompt_prefix1 = """You are a clinical geneticist. You have assembled known gene-disease associations from authoritative sources that include OMIM, GenCC, ClinGen, ClinVar, PanelApp, Orphanet, and dbNSFP. Now, you need to condense the phenotypes described in these different sources into a single concise comma-separate list that covers the primary features or symptoms of the disease, as well as the main organ systems that are affected. For example, when the source phenotypes are:

OMIM: 'Congenital disorder of glycosylation, type Ie', CLINGEN: 'congenital disorder of glycosylation type 1E', PANEL APP UK: 'Congenital disorder of glycosylation, type Ie, OMIM:608799, GDP-Man:Dol-P mannosyltransferase deficiency (Disorders of m
ultiple glycosylation and other glycosylation pathways); Congenital disorder of glycosylation, type Ie, OMIM:608799; Congenital disorder of glycosylation, type Ie, OMIM:608799; Congenital disorder of glycosylation, type Ie, OMIM:608799; Congenital
 disorder of glycosylation, type Ie, OMIM:608799, GDP-Man:Dol-P mannosyltransferase deficiency (Disorders of multiple glycosylation and other glycosylation pathways); Congenital disorder of glycosylation, type Ie, OMIM:608799, GDP-Man:Dol-P mannos
yltransferase deficiency (Disorders of multiple glycosylation and other glycosylation pathways); CONGENITAL DISORDERS OF GLYCOSYLATION; CONGENITAL DISORDERS OF GLYCOSYLATION 612379; Congenital disorder of glycosylation, type Ie, OMIM:608799; Conge
nital disorder of glycosylation, type Ie, OMIM:608799; Congenital disorder of glycosylation, type Ie, 608799', PANEL_APP_AU_phenotypes: 'Congenital disorder of glycosylation, type Ie, MIM# 608799; Congenital disorder of glycosylation, type Ie, 608
799; Congenital disorder of glycosylation, type Ie 608799; Congenital disorder of glycosylation, type Ie, 608799; Congenital disorder of glycosylation, type Ie, MIM# 608799; Congenital disorder of glycosylation, type Ie 608799; Congenital disorder
 of glycosylation, type Ie, 608799', DECIPHER: 'congenital disorder of glycosylation type 1E', CLINVAR: 'Congenital disorder of glycosylation type 1E, DPM1-related disorder, Inborn genetic diseases, not provided',
you would summarize this as:  "Congenital disorder of glycosylation type 1E".

The answer should contain only this summary and nothing else. There should be no intro or explanation - just the summary.
Now, try generating this type of summary for the following phenotype descriptions:

"""

# read the combined table. ncbi_gene_id is forced to str so its empty cells stay empty and its
# integer values aren't coerced to floats (e.g. "1.0") when written back out; load_bigquery.py then
# reads it as nullable Int64 so missing ids become BigQuery NULLs rather than 0.
df = pd.read_table(args.combined_table, dtype={"ncbi_gene_id": "str"})

# Phenotype columns in the combined table use "; " to join the per-gene, per-source values (see
# `separator` in generate_combined_gene_table.py). ClinGen/GenCC/PanelApp additionally carry a per-label
# classification/confidence column joined in the SAME order, so we can split both and keep only positives.
SEPARATOR = "; "

# IMPORTANT — POSITIVE ASSOCIATIONS ONLY: the LLM summary must describe phenotypes only from gene-disease
# associations that a source actually supports, never ones it disputes/refutes or rates as low-confidence.
# The combined table intentionally retains negative-evidence rows (their level is preserved in a parallel
# classification column), so this summarizer filters each source to its positive associations before
# prompting (drops ClinGen/GenCC/DECIPHER "Disputed"/"Refuted", non-"Green" PanelApp, OMIM Nondisease/
# Susceptibility, and Orphanet biomarker/susceptibility/fusion). OMIM provisional ("?") associations and
# Orphanet "candidate" associations are intentionally KEPT (still real, if tentative). This filtering is
# independent of what the pipeline loads into BigQuery / displays on the gene page for each source.
CLINGEN_POSITIVE_CLASSIFICATIONS = {"Definitive", "Strong", "Moderate", "Limited"}  # drop Disputed/Refuted/No Known Disease Relationship
GENCC_POSITIVE_CLASSIFICATIONS = {"Definitive", "Strong", "Moderate", "Limited", "Supportive"}  # drop Disputed/Refuted Evidence
DECIPHER_POSITIVE_CLASSIFICATIONS = {"Definitive", "Strong", "Moderate", "Limited", "Supportive"}  # DECIPHER GenCC validity; drop Disputed/Refuted
PANELAPP_POSITIVE_CONFIDENCE = {"3"}  # PanelApp "Green" (diagnostic-grade); drop 2=Amber, 1/0=Red
OMIM_DROP_CLASSIFICATIONS = {"Nondisease", "Susceptibility"}  # drop [] / {}; keep Confirmed, Provisional (?), and unclassified


def _has_content(v):
    """True if `v` carries an actual phenotype value.

    Treats NaN and separator/whitespace-only strings (e.g. "", "; ", ";;") as "no value". Genes that
    are on a PanelApp panel but have no OMIM phenotype produce an empty PANEL_APP_*_phenotypes value
    ("" or "; " after the per-source groupby join); pd.isna("") is False, so without this check such a
    gene would be treated as having a phenotype and trigger a pointless, billed LLM call on an empty
    prompt instead of falling back to the GWAS/constraint summary.
    """
    return not pd.isna(v) and str(v).strip(" ;,\t\r\n") != ""


def _positive_labels(labels_value, classifications_value, keep_classifications, strict=False):
    """Return the phenotype labels whose aligned classification is in `keep_classifications`.

    `labels_value` and `classifications_value` are the two "; "-joined columns for one gene, built by the
    same per-source groupby in generate_combined_gene_table.py, so token i of the labels pairs with token
    i of the classifications. On a token-count mismatch (e.g. a label that itself contains "; ", or a
    classification column missing entirely): with strict=False keep all non-empty labels (used where the
    labels are otherwise acceptable, e.g. OMIM); with strict=True drop everything (used where an unfiltered
    label could be a disputed/refuted association we must never emit, e.g. DECIPHER).
    """
    if not _has_content(labels_value):
        return []
    labels = [l.strip() for l in str(labels_value).split(SEPARATOR)]
    classifications = [c.strip() for c in str(classifications_value).split(SEPARATOR)]
    if len(labels) != len(classifications):
        return [] if strict else [l for l in labels if l]
    return [l for l, c in zip(labels, classifications) if l and c in keep_classifications]


def _labels_excluding(labels_value, classifications_value, drop_classifications):
    """Return the phenotype labels whose aligned classification is NOT in `drop_classifications`.

    Unlike _positive_labels (allowlist), this is a denylist: labels with an empty/unknown classification
    are KEPT. Used for OMIM, where we drop only Nondisease "[]" / Susceptibility "{}" and keep everything
    else (Confirmed, Provisional "?", and — on a cache written before OMIM_phenotype_classification existed
    — unclassified entries, so OMIM stays as-is until its cache is rebuilt).
    """
    if not _has_content(labels_value):
        return []
    labels = [l.strip() for l in str(labels_value).split(SEPARATOR)]
    classifications = [c.strip() for c in str(classifications_value).split(SEPARATOR)]
    if len(labels) != len(classifications):
        return [l for l in labels if l]
    return [l for l, c in zip(labels, classifications) if l and c not in drop_classifications]


def summarize_phenotypes(row):

    # Concatenate the phenotypes into a single string, by source — POSITIVE associations only (see the
    # module-level note). Each source is filtered by its own positive criterion so that what feeds the LLM
    # is independent of what the pipeline loads into BigQuery / displays on the gene page for that source:
    #   - OMIM: keep Confirmed + Provisional (drop Nondisease "[]" / Susceptibility "{}") via the aligned
    #     OMIM_phenotype_classification column.
    #   - ClinGen/GenCC/DECIPHER: keep positive validity classifications (drop Disputed/Refuted) via each
    #     source's aligned *_classification(s) column. DECIPHER uses strict=True so that if its
    #     classification column is ever missing/misaligned it contributes nothing rather than risk emitting
    #     a disputed/refuted association.
    #   - PanelApp UK/AU: keep only "Green" (confidence 3) via the aligned confidence column.
    #   - ClinVar (P/LP variants only), Fridman (curated), dbNSFP disease (UniProt involvement-in-disease)
    #     and dbNSFP HPO (HPO term NAMES, not ids) are inherently positive-only and pass through unfiltered.
    #   - Orphanet uses DBNSFP_orphanet_positive_disorder, pre-filtered upstream to genuine disease
    #     associations (candidate/modifier kept; biomarker/susceptibility/fusion dropped).
    phenotypes = []

    def add_source(label, value):
        if _has_content(value):
            phenotypes.append(f"{label}: {value}")

    add_source("OMIM", SEPARATOR.join(_labels_excluding(
        row.get("OMIM_phenotype_description"), row.get("OMIM_phenotype_classification"), OMIM_DROP_CLASSIFICATIONS)))
    # strict=True for the allowlist sources: their label columns retain disputed/refuted/non-Green
    # associations, so on any label↔classification token mismatch we must drop (never keep-all, which
    # could emit a negative association). Only OMIM (denylist above) keeps unclassified labels.
    add_source("CLINGEN", SEPARATOR.join(_positive_labels(
        row.get("CLINGEN_disease_label"), row.get("CLINGEN_classification"), CLINGEN_POSITIVE_CLASSIFICATIONS, strict=True)))
    add_source("GENCC", SEPARATOR.join(_positive_labels(
        row.get("GENCC_disease_name"), row.get("GENCC_classification"), GENCC_POSITIVE_CLASSIFICATIONS, strict=True)))
    add_source("DECIPHER", SEPARATOR.join(_positive_labels(
        row.get("DECIPHER_disease_names"), row.get("DECIPHER_classifications"), DECIPHER_POSITIVE_CLASSIFICATIONS, strict=True)))
    add_source("PANEL_APP_UK", SEPARATOR.join(_positive_labels(
        row.get("PANEL_APP_UK_phenotypes"), row.get("PANEL_APP_UK_confidence"), PANELAPP_POSITIVE_CONFIDENCE, strict=True)))
    add_source("PANEL_APP_AU", SEPARATOR.join(_positive_labels(
        row.get("PANEL_APP_AU_phenotypes"), row.get("PANEL_APP_AU_confidence"), PANELAPP_POSITIVE_CONFIDENCE, strict=True)))
    add_source("CLINVAR", row.get("CLINVAR_phenotypes"))
    add_source("FRIDMAN", row.get("FRIDMAN_phenotype_category"))
    add_source("ORPHANET", row.get("DBNSFP_orphanet_positive_disorder"))
    add_source("DBNSFP_DISEASE", row.get("DBNSFP_disease_description"))
    add_source("DBNSFP_HPO", row.get("DBNSFP_hpo_name"))

    if not phenotypes:
        if "GWAS_mondo_name" in row and _has_content(row["GWAS_mondo_name"]):
            return f"GWAS: " + str(row["GWAS_mondo_name"])
        else:
            constraint_type = []
            if row["pLI_v2"] >= 0.9:
                constraint_type.append("pLI_v2")
            if row["pLI_v4"] >= 0.9:
                constraint_type.append("pLI_v4")
            if row["lof_oe_ci_upper_v4"] <= 0.2:
                constraint_type.append("LOEUF")

            if constraint_type:
                return f"Constrained: {', '.join(constraint_type)}"
            else:
                return ""

    prompt = prompt_prefix1 + ", ".join(phenotypes)
    response = ask_gemini(prompt, model="2.5-flash", temperature=0, max_tokens=1000, system_prompt="", verbose=True)

    return response

# add a column for the phenotype summary
df["LLM_phenotype_summary"] = df.apply(summarize_phenotypes, axis=1)

# move the LLM_phenotype_summary column to be after the 'inheritance' column
initial_columns = [
    "ensembl_gene_id", "hgnc_gene_id", "refseq_id", "ncbi_gene_id", "gene_symbol", "gene_aliases",
    "in_MANE", "MANE_canonical_transcript_refseq_id", "MANE_canonical_transcript_ensembl_id",
    "MANE_clinical_transcript_refseq_id", "MANE_clinical_transcript_ensembl_id",
    "pLI_v2", "pLI_v4", "lof_oe_ci_upper_v4", "mis_oe_ci_upper_v4", "s_het",
    "inheritance",  "LLM_phenotype_summary", "sources",
    "gene_chrom", "gene_start", "gene_end",
]

df = df[initial_columns + [c for c in df.columns if c not in initial_columns]]

df.to_csv(args.output_path, sep="\t", index=False)

print(f"Wrote {len(df):,d} rows to {args.output_path}")
