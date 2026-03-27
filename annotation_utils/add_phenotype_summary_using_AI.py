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

# read the combined table
df = pd.read_table(args.combined_table)

def summarize_phenotypes(row):

    # concatenate the phenotypes into a single string, by source
    phenotypes = []
    for label, phenotype_column in [
        ("OMIM", "OMIM_phenotype_description"),
        ("CLINGEN", "CLINGEN_disease_label"),
        ("GENCC", "GENCC_disease_name"),
        ("PANEL_APP_UK", "PANEL_APP_UK_phenotypes"),
        ("PANEL_APP_AU", "PANEL_APP_AU_phenotypes"),
        ("CLINVAR", "CLINVAR_phenotypes"),
        ("FRIDMAN", "FRIDMAN_phenotype_category"),
        ("ORPHANET", "DBNSFP_orphanet_disorder"),
        ("DBNSFP_DISEASE", "DBNSFP_disease_description"),
    ]:
        if phenotype_column in row and not pd.isna(row[phenotype_column]):
            phenotypes.append(f"{label}: {row[phenotype_column]}")

    if not phenotypes:
        if "GWAS_mondo_name" in row and not pd.isna(row["GWAS_mondo_name"]):
            return f"GWAS: " + str(row["GWAS_mondo_name"])
        else:
            constraint_type = []
            if row["pLI_v2"] >= 0.9:
                constraint_type.append("pLI_v2")
            if row["pLI_v4"] >= 0.9:
                constraint_type.append("pLI_v4")
            if row["lof_oe_ci_upper_v4"] <= 0.2:
                constraint_type.append("LOEUF")
            if row["mis_oe_ci_upper_v4"] <= 0.2:
                constraint_type.append("MOEUF")

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
    "ensembl_gene_id", "hgnc_gene_id", "gene_symbol", "gene_aliases",  "pLI_v2", "pLI_v4", "lof_oe_ci_upper_v4", "mis_oe_ci_upper_v4", "s_het",
    "inheritance",  "LLM_phenotype_summary", "sources",
    "gene_chrom", "gene_start", "gene_end",
]

df = df[initial_columns + [c for c in df.columns if c not in initial_columns]]

df.to_csv(args.output_path, sep="\t", index=False)

print(f"Wrote {len(df):,d} rows to {args.output_path}")
