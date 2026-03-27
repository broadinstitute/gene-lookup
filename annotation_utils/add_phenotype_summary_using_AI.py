import argparse
import glob
import pandas as pd
from dotenv import load_dotenv
import time
import os
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

prompt_prefix2 = """
You are a clinical geneticist. You have assembled known gene-disease associations from authoritative sources that include OMIM, GenCC, ClinGen, ClinVar, PanelApp, Orphanet, and dbNSFP. Now, you need to select a single  
disease category that is the best match for the provided phenotypes. The possible disease categories are:

'BIOCHEMICAL/METABOLIC',
'CANCER',
'CARDIOVASCULAR',
'DEAFNESS',
'HEMATOLOGICAL',
'IMMUNOLOGICAL',
'OPHTHALMOLOGIC',
'NEUROLOGICAL',
'NEPHROLOGIC',
'PULMONARY',
'SKELETAL',
'RHEUMATOLOGIC & AUTOIMMUNE',
'ENDOCRINE',
'PRENATAL/REPRODUCTIVE',
'GASTROINTESTINAL',
'DYSMORPHOLOGY',
'SYNDROMIC',
'CONNECTIVE TISSUE',
'DERMATOLOGIC',
'PSYCHIATRIC'
If the phenotypes do not fit in one of these categories, put 'OTHER', or define a new disease area and label it.
In all scenarios, please output only the disease category and nothing else. There should be no intro or explanation - just the category.

For example, if the phenotypes are: 'Congenital disorder of glycosylation, type Ie, OMIM:608799, GDP-Man:Dol-P mannosyltransferase deficiency (Disorders of multiple glycosylation and other glycosylation pathways); Congenital disorder of glycosylation, type Ie, OMIM:608799, GDP-Man:Dol-P mannosyltransferase deficiency (Disorders of multiple glycosylation and other glycosylation pathways); CONGENITAL DISORDERS OF GLYCOSYLATION; CONGENITAL DISORDERS OF GLYCOSYLATION 612379; Congenital disorder of glycosylation, type Ie, OMIM:608799; Conge
nital disorder of glycosylation, type Ie, OMIM:608799; Congenital disorder of glycosylation, type Ie, 608799', you would output: 'BIOCHEMICAL/METABOLIC'.

or, if the phenotypes are: 'Microcephaly, short stature, and intellectual disability', you would output: 'DYSMORPHOLOGY'.

Now, please tell me the best disease category for the following phenotypes:
"""



# read the combined table
df = pd.read_table(args.combined_table)

def summarize_phenotypes(row, prompt_prefix=prompt_prefix1, blank_if_no_phenotypes=False):

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
        if blank_if_no_phenotypes:
            return ""
        
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

    prompt = prompt_prefix + ", ".join(phenotypes)
    response = ask_gemini(prompt, model="2.5-flash", temperature=0, max_tokens=1000, system_prompt="", verbose=True)

    return response

# add a column for the phenotype summary
df["LLM_phenotype_summary"] = df.apply(summarize_phenotypes, axis=1, prompt_prefix=prompt_prefix1, blank_if_no_phenotypes=False)
df["disease_category"] = df.apply(summarize_phenotypes, axis=1, prompt_prefix=prompt_prefix2, blank_if_no_phenotypes=True)


clingen_curations_dir = os.path.dirname(os.path.dirname(__file__))
clingen_curations_files = sorted(glob.glob(os.path.join(clingen_curations_dir, "CLINGEN_curations_export_at_*.csv")))
if not clingen_curations_files:
    print(f"WARNING: No CLINGEN_curations_export_at_*.csv files found in {clingen_curations_dir}. Skipping ClinGen curation processing.")
    df["clingen_curation"] = ""
else:
    clingen_curations_path = clingen_curations_files[-1]  # use the latest file
    print(f"Using ClinGen curations file: {clingen_curations_path}")
    clingen_df = pd.read_csv(clingen_curations_path)
    """Columns:
    Gene Symbol	Expert Panel	Curator	Disease Entity	Curation Type	Rationales	Uploaded Date	Precuration Date	Disease entity assigned Date	Precuration Complete Date	Curation Provisional Date	Curation Approved Date	Recuration assigned Date	Retired Assignment Date	Published Date	Classification	Created	GCI UUID
    """

    gene_name_to_clingen_curation_value = {}
    gene_alias_to_gene_name = {}
    duplicate_aliases = []
    for _, row in df.iterrows():
        if not pd.isna(row["gene_symbol"]) and row["gene_symbol"]:
            gene_name_to_clingen_curation_value[row["gene_symbol"].upper()] = None

        if not pd.isna(row["gene_aliases"]) and row["gene_aliases"] and not pd.isna(row["gene_symbol"]) and row["gene_symbol"]:
            gene_symbol = row["gene_symbol"].upper()
            for alias in row["gene_aliases"].split(","):
                alias = alias.strip().upper()
                if alias in gene_alias_to_gene_name and gene_alias_to_gene_name[alias] != gene_symbol:
                    #print(f"WARNING: Duplicate alias: {alias} for {row['gene_symbol']} and {gene_alias_to_gene_name[alias]}")
                    duplicate_aliases.append(alias)
                else:
                    gene_alias_to_gene_name[alias] = gene_symbol

    for alias in duplicate_aliases:
        if alias in gene_alias_to_gene_name:
            del gene_alias_to_gene_name[alias]

    clingen_gene_names_not_in_df = []
    for _, row in clingen_df.iterrows():
        if not pd.isna(row["Curation Type"]) and row["Curation Type"]:
            value = f"Curated: {row['Curation Type']}"
        else:
            value = "In Scope"

        gene_name = row["Gene Symbol"].strip().upper()
        gene_name = gene_alias_to_gene_name.get(gene_name, gene_name)
        if gene_name not in gene_name_to_clingen_curation_value:
            clingen_gene_names_not_in_df.append(gene_name)
        elif gene_name_to_clingen_curation_value[gene_name] is None or not gene_name_to_clingen_curation_value[gene_name].startswith("Curated"):
            gene_name_to_clingen_curation_value[gene_name] = value
        elif gene_name_to_clingen_curation_value[gene_name] != value and value is not None and value.startswith("Curated"):
            gene_name_to_clingen_curation_value[gene_name] += f", {value.replace('Curated: ', '')}"

    print(f"{len(clingen_gene_names_not_in_df):,} CLINGEN gene names not in df: {', '.join(clingen_gene_names_not_in_df)}")

    df["clingen_curation"] = df["gene_symbol"].str.upper().map(gene_name_to_clingen_curation_value)

# move the LLM_phenotype_summary column to be after the 'inheritance' column
initial_columns = [
    "ensembl_gene_id", "hgnc_gene_id", "gene_symbol", "gene_aliases",  "pLI_v2", "pLI_v4", "lof_oe_ci_upper_v4", "mis_oe_ci_upper_v4", "s_het",
    "inheritance",  "disease_category", "clingen_curation", "LLM_phenotype_summary", "sources",
    "gene_chrom", "gene_start", "gene_end",
]

df = df[initial_columns + [c for c in df.columns if c not in initial_columns]]

df.to_csv(args.output_path, sep="\t", index=False)

print(f"Wrote {len(df):,d} rows to {args.output_path}")

# print ClinGen stats
for label, df_subset in (
    ("Curated", df[df.clingen_curation.str.startswith("Curated") & df.clingen_curation.notna()]),
    ("In Scope", df[df.clingen_curation == "In Scope"]),
    ("Not in scope", df[df.clingen_curation.isna()]),
):
    df_subset = df_subset[df_subset.disease_category.notna() & (df_subset.disease_category.str.strip() != "")]

    print("-" * 100)
    print(f"{label}: {len(df_subset):,d} rows")
    for category, count in df_subset.disease_category.value_counts().items():
        print(f"{count:,d}\t{category}")
    print()
