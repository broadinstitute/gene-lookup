import argparse
from datetime import datetime
import os
import pandas as pd

from annotation_utils.get_omim_table import get_omim_table
from annotation_utils.get_clingen_table import get_clingen_gene_disease_validity_table, get_clingen_haploinsufficient_genes_table
from annotation_utils.get_hgnc_table import get_hgnc_table, get_hgnc_to_ensg_id_map, get_ensg_id_to_hgnc_id_map
from annotation_utils.get_panel_app_table import get_panel_app_table
from annotation_utils.get_ensembl_db_info import get_transcript_id_to_gene_id, get_gene_metadata
from annotation_utils.get_gwas_catalog import get_gwas_catalog_rare_disease_records
from annotation_utils.get_gencc_table import get_gencc_table
from annotation_utils.get_decipher_genes import get_decipher_gene_table
from annotation_utils.get_clinvar_table import get_clinvar_gene_disease_table
from annotation_utils.get_constraint_scores import get_constraint_scores

parser = argparse.ArgumentParser()
parser.add_argument("--skip-gwas", action="store_true", help="Don't add columns related to GWAS catalog rare disease records")
parser.add_argument("--skip-fridman", action="store_true", help="Don't add columns related to the Fridman et al. 2025 list of recessive genes")
parser.add_argument("--print-example-genes", action="store_true", help="Print example genes")
parser.add_argument("--force", action="store_true", help="Force re-download of all data, even if cached")
args = parser.parse_args()

if args.force:
    os.environ["FORCE_DOWNLOAD"] = "1"

include_GWAS = not args.skip_gwas
include_Fridman = not args.skip_fridman
print_example_genes = args.print_example_genes

separator = "; "

def safe_max(values, default=float('nan')):
    """Return max of non-empty values, or default if all values are empty."""
    non_empty = [float(v) for v in values if v != ""]
    return max(non_empty) if non_empty else default

def safe_min(values, default=float('nan')):
    """Return min of non-empty values, or default if all values are empty."""
    non_empty = [float(v) for v in values if v != ""]
    return min(non_empty) if non_empty else default

def normalize_nulls(x):
    x = x if not pd.isna(x) else ""
    if isinstance(x, float) and x.is_integer():
        x = int(x)

    return str(x)

HGNC_to_ENSG_map = get_hgnc_to_ensg_id_map()
ENSG_to_HGNC_map = get_ensg_id_to_hgnc_id_map()
df_hgnc = get_hgnc_table()
ENSG_to_gene_name_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["Approved symbol"]))
ENSG_to_gene_name_aliases_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["Alias symbols"]))

"""
Example:

chrom                                                   1
start                                             1013497
end                                               1014540
mim_number                                         147571
phenotype_mim_number                               616126
phenotypic_series_number
phenotype_inheritance                 Autosomal recessive
gene_symbols                  G1P2,  IFI15,  IMD38, ISG15
gene_id                                   ENSG00000187608
gene_description            ISG15 ubiquitin-like modifier
phenotype_description                 Immunodeficiency 38
date_created
date_updated
mouse_gene_id                         Isg15 (MGI:1855694)
oe_lof_upper                                        1.691
pLI                                               0.40527
mis_z                                           -0.044129
"""

print("Getting OMIM table")
df_omim = get_omim_table()
print(f"Got {len(df_omim):,d} rows from OMIM")
df_omim = df_omim[df_omim["phenotype_mim_number"].notna() & (df_omim["phenotype_mim_number"] != "")]
print(" "*8, f"Kept {len(df_omim):,d} rows from OMIM, containing {len(df_omim['gene_id'].unique()):,d} unique genes")
df_omim = df_omim[[
    "gene_id",   # ENSG id
    "mim_number",
    "phenotype_mim_number",
    "phenotype_inheritance",
    "phenotype_description",
    #"phenotypic_series_number",
    #"gene_symbols",
    #"gene_description",
    #"mouse_gene_id",
    #"oe_lof_upper",
    #"pLI",
    #"mis_z",
]]

df_omim = df_omim.rename(columns={
    "gene_id": "OMIM_gene_id",
    "mim_number": "OMIM_mim_number",
    "phenotype_mim_number": "OMIM_phenotype_mim_number",
    "phenotype_inheritance": "OMIM_inheritance",
    "phenotype_description": "OMIM_phenotype_description",
    #"phenotypic_series_number": "OMIM_phenotypic_series_number",
    #"oe_lof_upper": "LOEUF",
    #"pLI": "pLI",
    #"mis_z": "mis_z",
})


before = len(df_omim)
df_omim = df_omim[df_omim["OMIM_gene_id"].notna() & (df_omim["OMIM_phenotype_mim_number"].notna() | df_omim["OMIM_phenotype_description"].notna())]
df_omim = df_omim[(df_omim["OMIM_gene_id"] != "") & ((df_omim["OMIM_phenotype_mim_number"] != "") | (df_omim["OMIM_phenotype_description"] != ""))]
print("\t", f"Kept {len(df_omim):,d} out of {before:,d} ({(len(df_omim) / before):.1%}) rows which had both a gene and a phenotype")


# group by gene_id and combine the other fields using ; as a separator
df_omim = df_omim.groupby("OMIM_gene_id").agg({
    "OMIM_mim_number": lambda x: separator.join(normalize_nulls(v) for v in x),
    "OMIM_phenotype_mim_number": lambda x: separator.join(normalize_nulls(v) for v in x),
    "OMIM_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
    "OMIM_phenotype_description": lambda x: separator.join(normalize_nulls(v) for v in x),
    #"OMIM_phenotypic_series_number": lambda x: separator.join(normalize_nulls(v) for v in x),
    #"LOEUF": lambda x: normalize_nulls(x.iloc[0]),
    #"pLI": lambda x: normalize_nulls(x.iloc[0]),
    #"mis_z": lambda x: normalize_nulls(x.iloc[0]),
}).reset_index()


df_omim.set_index("OMIM_gene_id", inplace=True)



print("Getting ClinGen haploinsufficient genes table")
df_clingen_haploinsufficient_genes = get_clingen_haploinsufficient_genes_table()
df_clingen_haploinsufficient_genes = df_clingen_haploinsufficient_genes.rename(columns={
    "HAPLOINSUFFICIENCY": "CLINGEN_haploinsufficient",
})

df_clingen_haploinsufficient_genes = df_clingen_haploinsufficient_genes[[
    "HGNC ID",
    "CLINGEN_haploinsufficient",
]]

df_clingen_haploinsufficient_genes.set_index("HGNC ID", inplace=True)


"""
Example:

GENE SYMBOL                                                        AARS1
GENE ID (HGNC)                                                   HGNC:20
DISEASE LABEL                 Charcot-Marie-Tooth disease axonal type 2N
DISEASE ID (MONDO)                                         MONDO:0013212
MOI                                                                   AD
SOP                                                                SOP10
CLASSIFICATION                                                Definitive
ONLINE REPORT          https://search.clinicalgenome.org/kb/gene-vali...
CLASSIFICATION DATE                             2024-03-14T16:00:00.000Z
GCEP                   Charcot-Marie-Tooth Disease Gene Curation Expe...
"""
print("Getting ClinGen gene disease validity table")
df_clingen = get_clingen_gene_disease_validity_table()
print(f"Got {len(df_clingen):,d} rows from ClinGen, containing {len(df_clingen['GENE ID (HGNC)'].unique()):,d} unique genes")


df_clingen = df_clingen.rename(columns={
    "DISEASE LABEL": "CLINGEN_disease_label",
    "DISEASE ID (MONDO)": "CLINGEN_disease_mondo_id",
    "MOI": "CLINGEN_inheritance",
    "CLASSIFICATION": "CLINGEN_classification",
})

df_clingen = df_clingen.set_index("GENE ID (HGNC)").join(df_clingen_haploinsufficient_genes, how="outer").reset_index()

# Convert HGNC -> ENSG for joining, but preserve native HGNC ID
df_clingen["CLINGEN_gene_id"] = df_clingen["GENE ID (HGNC)"].map(HGNC_to_ENSG_map)
df_clingen["CLINGEN_hgnc_gene_id"] = df_clingen["GENE ID (HGNC)"]
hgnc_ids_with_missing_esng = df_clingen[df_clingen['CLINGEN_gene_id'].isna()]['GENE ID (HGNC)'].unique()
assert len(hgnc_ids_with_missing_esng) == 0, f"Could not convert the following HGNC ids to ENSG: {', '.join(hgnc_ids_with_missing_esng)}"

df_clingen = df_clingen[[
    "CLINGEN_gene_id",
    "CLINGEN_hgnc_gene_id",
    "CLINGEN_disease_label",
    "CLINGEN_disease_mondo_id",
    "CLINGEN_inheritance",
    "CLINGEN_classification",
    "CLINGEN_haploinsufficient",
]]


before = len(df_clingen)
df_clingen = df_clingen[df_clingen["CLINGEN_classification"].isin({
    "Definitive", "Limited", "Moderate", "Strong"
})]
print("\t", f"Kept {len(df_clingen):,d} out of {before:,d} ({(len(df_clingen) / before):.1%}) rows which had a Definitive, Limited, Moderate, or Strong classification")



# group by CLINGEN_gene_id and combine the other fields using ; as a separator
df_clingen = df_clingen.groupby("CLINGEN_gene_id").agg({
    "CLINGEN_hgnc_gene_id": "first",
    "CLINGEN_disease_label": lambda x: separator.join(normalize_nulls(v) for v in x),
    "CLINGEN_disease_mondo_id": lambda x: separator.join(normalize_nulls(v) for v in x),
    "CLINGEN_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
    "CLINGEN_classification": lambda x: separator.join(normalize_nulls(v) for v in x),
    "CLINGEN_haploinsufficient": lambda x: separator.join(sorted(set(normalize_nulls(v) for v in x))),
}).reset_index()

df_clingen.set_index("CLINGEN_gene_id", inplace=True)


print("Getting constraint scores table")
df_constraint_scores = get_constraint_scores()

"""Example:

wm8c1-2cf:~/code/annotation-utils $ python3 annotation_utils/get_constraint_scores.py
                       pLI_v2    pLI_v4  lof_oe_v4  lof_oe_ci_lower_v4  lof_oe_ci_upper_v4  mis_oe_v4  mis_oe_ci_lower_v4  mis_oe_ci_upper_v4
gene_id
ENSG00000000003  6.607900e-02       NaN        NaN                 NaN                 NaN        NaN                 NaN                 NaN
ENSG00000000005  3.511900e-03       NaN        NaN                 NaN                 NaN        NaN                 NaN                 NaN
ENSG00000000419  1.623700e-04  0.124840    0.48616               0.265               0.960     0.8549               0.719               1.020
ENSG00000000457  2.815100e-01  0.549410    0.38525               0.250               0.613     0.8228               0.764               0.886
ENSG00000000460  3.352300e-10  0.003679    0.78335               0.459               1.408     1.0024               0.865               1.164
"""

df_constraint_scores = df_constraint_scores.rename(columns={
    "gene_id": "CONSTRAINT_scores_gene_id",
})

df_constraint_scores = df_constraint_scores[[
    "CONSTRAINT_scores_gene_id",
    "pLI_v2",
    "pLI_v4",
    #"lof_oe_v4",
    #"lof_oe_ci_lower_v4",
    "lof_oe_ci_upper_v4",
    #"mis_oe_v4",
    #"mis_oe_ci_lower_v4",
    "mis_oe_ci_upper_v4",
]]

df_constraint_scores = df_constraint_scores.groupby("CONSTRAINT_scores_gene_id").agg({
    "pLI_v2": lambda x: safe_max(x),
    "pLI_v4": lambda x: safe_max(x),
    "lof_oe_ci_upper_v4": lambda x: safe_min(x),
    "mis_oe_ci_upper_v4": lambda x: safe_min(x),
}).reset_index()

df_constraint_scores.set_index("CONSTRAINT_scores_gene_id", inplace=True)




"""
source                                               PanelApp UK
hgnc                                                  HGNC:24641
gene_name                    chromosome 16 open reading frame 62
biotype                                           protein_coding
gene_id                                          ENSG00000103544
confidence                                                     2
penetrance                                                  None
mode_of_pathogenicity                                       None
mode_of_inheritance      BIALLELIC, autosomal or pseudoautosomal
publications                                            31712251
evidence                         Expert Review Amber, Literature
phenotypes                    3C/Ritscher-Schinzel-like syndrome
panel_name                             Chondrodysplasia punctata
"""

print("Getting PanelApp tables")
df_panel_app = get_panel_app_table()
PANEL_APP_UK_LABEL = "PanelApp UK"
PANEL_APP_AU_LABEL = "PanelApp Australia"
assert set(df_panel_app["source"]) == {
    PANEL_APP_UK_LABEL, PANEL_APP_AU_LABEL,
}
print(f"Got {len(df_panel_app):,d} rows from PanelApp, containing {len(df_panel_app['gene_id'].unique()):,d} unique genes")

df_panel_app["phenotypes"] = df_panel_app["phenotypes"].str.replace("No OMIM phenotype", "")

before = len(df_panel_app)
df_panel_app["gene_id"] = df_panel_app["gene_id"].fillna(df_panel_app["hgnc"].map(HGNC_to_ENSG_map))
df_panel_app = df_panel_app[df_panel_app["gene_id"].notna() & (df_panel_app["gene_id"] != "")]
df_panel_app = df_panel_app[df_panel_app["phenotypes"].notna() & (df_panel_app["phenotypes"] != "")]
print("\t", f"Kept {len(df_panel_app):,d} out of {before:,d} ({(len(df_panel_app) / before):.1%}) rows which had a gene id and a phenotype")

before = len(df_panel_app)
df_panel_app = df_panel_app[~df_panel_app["evidence"].apply(lambda x: not isinstance(x, str) or "Expert Review Red" in x)]
print("\t", f"Kept {len(df_panel_app):,d} out of {before:,d} ({(len(df_panel_app) / before):.1%}) rows which had evidence other than 'Expert Review Red'")

df_panel_app = df_panel_app[[
    "gene_id",
    "source",
    #"hgnc",
    #"gene_name",
    #"biotype",
    "confidence",
    "penetrance",
    "mode_of_pathogenicity",
    "mode_of_inheritance",
    #"publications",
    "evidence",
    "phenotypes",
    "panel_name",
]]

df_panel_app = df_panel_app.groupby(["gene_id", "source"]).agg({
    "confidence": lambda x: separator.join(normalize_nulls(v) for v in x),
    "penetrance": lambda x: separator.join(normalize_nulls(v) for v in x),
    "mode_of_pathogenicity": lambda x: separator.join(normalize_nulls(v) for v in x),
    "mode_of_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
    "evidence": lambda x: separator.join(normalize_nulls(v) for v in x),
    "phenotypes": lambda x: separator.join(normalize_nulls(v) for v in x),
    "panel_name": lambda x: separator.join(normalize_nulls(v) for v in x),
}).reset_index()

# drop column "source"
df_panel_app_uk = df_panel_app[df_panel_app["source"] == PANEL_APP_UK_LABEL].drop("source", axis=1)
df_panel_app_uk["InPanelAppUK"] = True
df_panel_app_au = df_panel_app[df_panel_app["source"] == PANEL_APP_AU_LABEL].drop("source", axis=1)
df_panel_app_au["InPanelAppAU"] = True

for panel_app_label, df_pannel_app in [
    ("UK", df_panel_app_uk),
    ("AU", df_panel_app_au),
]:
    df_pannel_app.rename(columns={
        "gene_id": f"PANEL_APP_{panel_app_label}_gene_id",
        "confidence": f"PANEL_APP_{panel_app_label}_confidence",
        "penetrance": f"PANEL_APP_{panel_app_label}_penetrance",
        "mode_of_pathogenicity": f"PANEL_APP_{panel_app_label}_mode_of_pathogenicity",
        "mode_of_inheritance": f"PANEL_APP_{panel_app_label}_inheritance",
        "evidence": f"PANEL_APP_{panel_app_label}_evidence",
        "phenotypes": f"PANEL_APP_{panel_app_label}_phenotypes",
        "panel_name": f"PANEL_APP_{panel_app_label}_panel_name",
    }, inplace=True)

    df_pannel_app.set_index(f"PANEL_APP_{panel_app_label}_gene_id", inplace=True)

    print("\t", f"PanelApp {panel_app_label} contains {len(df_pannel_app):,d} gene ids")
    # group by gene_id and combine the other fields using ; as a separator

# do an outer join of the 2 tables
df_panel_app = pd.merge(df_panel_app_uk, df_panel_app_au, how="outer", left_index=True, right_index=True)
print("\t", f"Merged PanelApp table contains {len(df_panel_app):,d} gene ids")


print("Getting table of recessive genes from Fridman et al. 2025")
transcript_id_to_gene_id = get_transcript_id_to_gene_id()

fridman_path = "annotation_utils/data/AR_genes_from_Fridman_2025.tsv"
if include_Fridman and not os.path.exists(fridman_path):
    print(f"WARNING: Fridman file not found at {fridman_path}. Skipping Fridman data.")
    include_Fridman = False

if include_Fridman:
    df_fridman = pd.read_table(fridman_path)
    df_fridman["FRIDMAN_gene_id"] = df_fridman["Transcripts"].apply(
        lambda transcript_list: ", ".join(transcript_id_to_gene_id[t] for t in transcript_list.split(",") if t in transcript_id_to_gene_id)
    )
    assert sum(df_fridman["FRIDMAN_gene_id"].str.contains(",")) == 0, "Some rows had multiple gene ids: " + str(df_fridman[df_fridman["FRIDMAN_gene_id"].str.contains(",")])

    df_fridman = df_fridman[df_fridman["FRIDMAN_gene_id"].notna() & (df_fridman["FRIDMAN_gene_id"] != "")]
    df_fridman = df_fridman[[
        "FRIDMAN_gene_id",
        #"Transcripts",
        "OMIM phenotype ID",
        "Disorder group",
        "Inheritance mode (AR/AR-AD)"
    ]]

    df_fridman.rename(columns={
        "OMIM phenotype ID": "FRIDMAN_omim_phenotype_id",
        "Disorder group": "FRIDMAN_phenotype_category",
        "Inheritance mode (AR/AR-AD)": "FRIDMAN_inheritance",
    }, inplace=True)


    df_fridman = df_fridman.groupby("FRIDMAN_gene_id").agg({
        "FRIDMAN_omim_phenotype_id": lambda x: separator.join(normalize_nulls(v) for v in x),
        "FRIDMAN_phenotype_category": lambda x: separator.join(normalize_nulls(v) for v in x),
        "FRIDMAN_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
    }).reset_index()

    df_fridman["InFridman"] = True

    df_fridman.set_index("FRIDMAN_gene_id", inplace=True)


df_gwas = None
if include_GWAS:
    print("Getting GWAS catalog rare disease records")
    df_gwas = get_gwas_catalog_rare_disease_records()

    """
    MONDO_ID           MONDO:0016158
    CHR_ID                         2
    CHR_POS                241837710
    SNPS                  rs34071003
    P-VALUE                 0.000004
    OR or BETA                 1.328
    95% CI (TEXT)    [1.1789-1.4980]
    GENE_ID          ENSG00000204099
    GENE_TYPE               UPSTREAM
    GENE_DISTANCE            20297.0
    """

    df_gwas = df_gwas[~df_gwas["MONDO_CATEGORY"].isin({"cancer or benign tumor", "infectious disease"})]

    df_gwas.rename(columns={
        "MONDO_ID": "GWAS_mondo_id",
        "MONDO_NAME": "GWAS_mondo_name",
        "MONDO_CATEGORY": "GWAS_mondo_category",
        "CHR_ID": "GWAS_chr_id",
        "CHR_POS": "GWAS_chr_pos",
        "SNPS": "GWAS_snps",
        "P-VALUE": "GWAS_p_value",
        "OR or BETA": "GWAS_odds_ratio_or_beta",
        "95% CI (TEXT)": "GWAS_95_ci_text",
        "GENE_ID": "GWAS_gene_id",
        "GENE_TYPE": "GWAS_gene_type",
        "GENE_DISTANCE": "GWAS_gene_distance",
    }, inplace=True)

    df_gwas = df_gwas[df_gwas["GWAS_gene_id"].notna() & (df_gwas["GWAS_gene_id"] != "")]

    df_gwas = df_gwas.groupby("GWAS_gene_id").agg({
        "GWAS_mondo_id": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_mondo_name": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_mondo_category": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_gene_type": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_gene_distance": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_chr_id": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_chr_pos": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_snps": lambda x: separator.join(normalize_nulls(v) for v in x),
        "GWAS_p_value": lambda x: safe_min(x),
        "GWAS_odds_ratio_or_beta": lambda x: safe_min(x),
        "GWAS_95_ci_text": lambda x: separator.join(normalize_nulls(v) for v in x),
    }).reset_index()

    df_gwas["InGWAS"] = True
    df_gwas.set_index("GWAS_gene_id", inplace=True)


print("Getting GenCC table")
df_gencc = get_gencc_table()

"""
[
    "gene_id",
    "hgnc_gene_id",
    "disease_id",
    "disease_name",
    "classification",
    "inheritance",
]
"""
df_gencc.rename(columns={
    "gene_id": "GENCC_gene_id",
    "hgnc_gene_id": "GENCC_hgnc_gene_id",
    "disease_id": "GENCC_disease_id",
    "disease_name": "GENCC_disease_name",
    "classification": "GENCC_classification",
    "inheritance": "GENCC_inheritance",
}, inplace=True)

df_gencc = df_gencc[[
    "GENCC_gene_id",
    "GENCC_hgnc_gene_id",
    #"GENCC_disease_id",
    "GENCC_disease_name",
    "GENCC_classification",
    "GENCC_inheritance",
]]

df_gencc = df_gencc.groupby("GENCC_gene_id").agg({
    "GENCC_hgnc_gene_id": "first",
    "GENCC_disease_name": lambda x: separator.join(normalize_nulls(v) for v in x),
    "GENCC_classification": lambda x: separator.join(normalize_nulls(v) for v in x),
    "GENCC_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
}).reset_index()

df_gencc.set_index("GENCC_gene_id", inplace=True)

print("Getting Decipher table")
df_decipher = get_decipher_gene_table()

df_decipher.rename(columns={
    "gene_id": "DECIPHER_gene_id",
    "inheritance_modes": "DECIPHER_inheritance",
    "disease_names": "DECIPHER_disease_names",
}, inplace=True)

df_decipher = df_decipher.groupby("DECIPHER_gene_id").agg({
    "DECIPHER_inheritance": lambda x: separator.join(normalize_nulls(v) for v in x),
    "DECIPHER_disease_names": lambda x: separator.join(normalize_nulls(v) for v in x),
}).reset_index()

df_decipher.set_index("DECIPHER_gene_id", inplace=True)

print("Getting ClinVar table")
df_clinvar = get_clinvar_gene_disease_table()

df_clinvar.rename(columns={
    "gene_id": "CLINVAR_gene_id",
    "phenotypes": "CLINVAR_phenotypes",
    "clinical_significance": "CLINVAR_clinical_significance",
    "gold_stars": "CLINVAR_stars",
    "major_consequences": "CLINVAR_variant_consequences",
}, inplace=True)

df_clinvar.set_index("CLINVAR_gene_id", inplace=True)


print(f"Merging "
      f"OMIM ({len(df_omim):,d} rows) "
      f"ClinGen ({len(df_clingen):,d} rows) "
      f"PanelApp ({len(df_panel_app):,d} rows) "
      f"GenCC ({len(df_gencc):,d} rows) "
      f"Decipher ({len(df_decipher):,d} rows) "
      f"ClinVar ({len(df_clinvar):,d} rows) " +
      (f"GWAS catalog ({len(df_gwas):,d} rows) " if include_GWAS else "") +
      (f"Fridman ({len(df_fridman):,d} rows) " if include_Fridman else "")
)

df_omim["InOMIM"] = True
df_clingen["InClinGen"] = True
df_gencc["InGenCC"] = True
df_decipher["InDecipher"] = True
df_clinvar["InClinVar"] = True


before = list(df_omim.index)
print(f"Starting with {len(df_omim):,d} genes from OMIM")
df_combined = pd.merge(df_omim, df_clingen, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging OMIM and ClinGen"
print(f"Added {len(df_combined) - len(before):,d} genes from ClinGen" + (f" - examples: {', '.join(list(sorted(set(df_combined.index) - set(before)))[:5])}" if print_example_genes else ""))

before = list(df_combined.index)
df_combined = pd.merge(df_combined, df_gencc, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with GenCC"
print(f"Added {len(df_combined) - len(before):,d} genes from GenCC" + (f" - examples: {', '.join(list(sorted(set(df_combined.index) - set(before)))[:5])}" if print_example_genes else ""))

before = list(df_combined.index)
df_combined = pd.merge(df_combined, df_panel_app, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with PanelApp"
print(f"Added {len(df_combined) - len(before):,d} genes from PanelApp" + (f" - examples: {', '.join(list(sorted(set(df_combined.index) - set(before)))[:5])}" if print_example_genes else ""))

before = list(df_combined.index)
df_combined = pd.merge(df_combined, df_decipher, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with Decipher"
print(f"Added {len(df_combined) - len(before):,d} genes from Decipher" + (f" - examples: {', '.join(list(sorted(set(df_combined.index) - set(before)))[:5])}" if print_example_genes else ""))

before = list(df_combined.index)
df_combined = pd.merge(df_combined, df_clinvar, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with ClinVar"
print(f"Added {len(df_combined) - len(before):,d} genes from ClinVar" + (f" - examples: {', '.join(list(set(sorted(df_combined.index)) - set(sorted(before)))[:20])}" if print_example_genes else ""))

if include_GWAS:
    before = list(df_combined.index)
    df_combined = pd.merge(df_combined, df_gwas, how="outer", left_index=True, right_index=True)
    assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with GWAS"
    print(f"Added {len(df_combined) - len(before):,d} genes from GWAS catalog" + (f" - examples: {', '.join(list(sorted(set(df_combined.index) - set(before)))[:5])}" if print_example_genes else ""))

if include_Fridman:
    before = list(df_combined.index)
    df_combined = pd.merge(df_combined, df_fridman, how="outer", left_index=True, right_index=True)
    assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with Fridman"
    print(f"Added {len(df_combined) - len(before):,d} genes from Fridman et al. 2025 list of recessive genes" + (f" - examples: {', '.join(list(set(df_combined.index) - set(before))[:5])}" if print_example_genes else ""))

# add gene chrom, start, end
df_gene_chrom_start_end = pd.DataFrame(get_gene_metadata().values())
df_gene_chrom_start_end = df_gene_chrom_start_end[["gene.stable_id", "chrom", "start", "end"]]
df_gene_chrom_start_end.set_index("gene.stable_id", inplace=True)
missing_gene_ids = set(df_combined.index) - set(df_gene_chrom_start_end.index)
if len(missing_gene_ids) > 0:
    print(f"WARNING: chrom/start/end not available for {len(missing_gene_ids):,d} genes: {', '.join(list(missing_gene_ids)[:20])}" + (", ..." if len(missing_gene_ids) > 20 else ""))
df_combined = pd.merge(df_combined, df_gene_chrom_start_end, how="left", left_index=True, right_index=True)
df_combined["start"] = df_combined["start"].fillna(0).astype(int)
df_combined["end"] = df_combined["end"].fillna(0).astype(int)

# add constraint scores
before = list(df_combined.index)
df_combined = pd.merge(df_combined, df_constraint_scores, how="left", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with constraint scores"
rows_with_constraint_scores = sum(df_combined["pLI_v2"].notna() | df_combined["pLI_v4"].notna() | df_combined["lof_oe_ci_upper_v4"].notna() | df_combined["mis_oe_ci_upper_v4"].notna())
print(f"Added constraint scores to {rows_with_constraint_scores:,d} out of {len(df_combined):,d} ({(rows_with_constraint_scores / len(df_combined)):.1%}) genes")

pLI_v2_THRESHOLD = 0.9
pLI_v4_THRESHOLD = 0.9
LOEUF_CONSTRAINT_THRESHOLD = 0.2
MOEUF_CONSTRAINT_THRESHOLD = 0.2

# add highly constrained genes that are not in the other sources:
df_highly_constrained_genes = df_constraint_scores[
    ~df_constraint_scores.index.isin(df_combined.index) & (
        (df_constraint_scores["pLI_v2"] >= pLI_v2_THRESHOLD) |
        (df_constraint_scores["pLI_v4"] >= pLI_v4_THRESHOLD) |
        (df_constraint_scores["lof_oe_ci_upper_v4"] <= LOEUF_CONSTRAINT_THRESHOLD) |
        (df_constraint_scores["mis_oe_ci_upper_v4"] <= MOEUF_CONSTRAINT_THRESHOLD)
    )
]

print(f"Added {len(df_highly_constrained_genes):,d} highly constrained genes that are not in the other sources")
df_combined = pd.concat([df_combined, df_highly_constrained_genes])


df_combined.reset_index(inplace=True)
df_combined.rename(columns={
    "index": "ensembl_gene_id",
}, inplace=True)

# Build hgnc_gene_id by coalescing: native HGNC from ClinGen/GenCC first, then ENSG->HGNC map as fallback
df_combined["hgnc_gene_id"] = df_combined["CLINGEN_hgnc_gene_id"].combine_first(
    df_combined["GENCC_hgnc_gene_id"]
).combine_first(
    df_combined["ensembl_gene_id"].map(ENSG_to_HGNC_map)
)
df_combined.drop(columns=["CLINGEN_hgnc_gene_id", "GENCC_hgnc_gene_id"], inplace=True)

df_combined["gene_symbol"] = df_combined["ensembl_gene_id"].map(ENSG_to_gene_name_map).str.upper()
df_combined["gene_aliases"] = df_combined["ensembl_gene_id"].map(ENSG_to_gene_name_aliases_map).str.upper()
#if df_combined["hgnc_gene_id"].isna().sum() > 0:
#    print(f"WARNING: {df_combined['hgnc_gene_id'].isna().sum():,d} genes had no HGNC id")
#    print(df_combined[df_combined["hgnc_gene_id"].isna()])

def compute_sources_string(row):
    sources = []
    if row["InOMIM"] == True:
        sources.append("OMIM")
    if row["InClinGen"] == True:
        sources.append("ClinGen")
    if row["InGenCC"] == True:
        sources.append("GenCC")
    if row["InPanelAppUK"] == True:
        sources.append("PanelAppUK")
    if row["InPanelAppAU"] == True:
        sources.append("PanelAppAU")
    if row["InDecipher"] == True:
        sources.append("Decipher")
    if row["InClinVar"] == True:
        sources.append("ClinVar")
    if include_GWAS and row["InGWAS"] == True:
        sources.append("GWAS")
    if include_Fridman and row["InFridman"] == True:
        sources.append("Fridman")
    if row["pLI_v2"] >= pLI_v2_THRESHOLD:
        sources.append("pLI_v2:" + str(round(row["pLI_v2"], 2)))
    if row["pLI_v4"] >= pLI_v4_THRESHOLD:
        sources.append("pLI_v4:" + str(round(row["pLI_v4"], 2)))
    if row["lof_oe_ci_upper_v4"] <= LOEUF_CONSTRAINT_THRESHOLD:
        sources.append("LOEUF:" + str(round(row["lof_oe_ci_upper_v4"], 2)))
    if row["mis_oe_ci_upper_v4"] <= MOEUF_CONSTRAINT_THRESHOLD:
        sources.append("MOEUF:" + str(round(row["mis_oe_ci_upper_v4"], 2)))

    return f"{len(sources)}: " + ", ".join(sources)

df_combined["sources"] = df_combined.apply(compute_sources_string, axis=1)
df_combined.drop(columns=["InOMIM", "InClinGen", "InGenCC", "InPanelAppAU", "InPanelAppUK", "InDecipher", "InClinVar"], inplace=True)

if include_GWAS:
    df_combined.drop(columns=["InGWAS"], inplace=True)
if include_Fridman:
    df_combined.drop(columns=["InFridman"], inplace=True)

def summarize_inheritance(row):
    inheritance = set()
    columns = [
        "OMIM_inheritance",
        "CLINGEN_inheritance",
        "GENCC_inheritance",
        "PANEL_APP_UK_inheritance",
        "PANEL_APP_AU_inheritance",
        "DECIPHER_inheritance",
    ]
    if include_Fridman:
        columns.append("FRIDMAN_inheritance")
    for column in columns:
        if column in row and isinstance(row[column], str):
            inheritance.update({i.strip() for i in row[column].split(";") if i.strip() != ""})

    return "" if len(inheritance) == 0 else "; ".join(sorted(inheritance))

df_combined["inheritance"] = df_combined.apply(summarize_inheritance, axis=1)

# move the gene_id, hgnc_gene_id, gene_symbol, and gene_aliases columns to the front
initial_columns = ["ensembl_gene_id", "hgnc_gene_id", "gene_symbol", "gene_aliases", "pLI_v2", "pLI_v4", "lof_oe_ci_upper_v4", "mis_oe_ci_upper_v4", "inheritance", "sources"]
df_combined = df_combined[initial_columns + [c for c in df_combined.columns if c not in initial_columns]]
#df_combined.sort_values(by=["sources", "gene_id"], inplace=True)

#timestamp = datetime.now().strftime("%Y_%m_%d")
#output_path = f"combined_mendelian_gene_disease_table.{len(df_combined)}_genes.{timestamp}.tsv"
output_path = f"combined_mendelian_gene_disease_table.tsv.gz"
df_combined.to_csv(output_path, sep="\t", index=False)
print(f"Wrote {len(df_combined):,d} genes to {output_path}")

output_path = f"combined_mendelian_gene_disease_table.only_in_clinvar.tsv.gz"
df_clinvar_only = df_combined[(df_combined["sources"] == "1: ClinVar") | (df_combined["sources"] == "2: ClinVar, Fridman")]
df_clinvar_only.to_csv(output_path, sep="\t", index=False)
print(f"Wrote {len(df_clinvar_only):,d} genes to {output_path}")
