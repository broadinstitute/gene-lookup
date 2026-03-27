import io
import os

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table

DBNSFP_BASE_URL = "https://dist.genos.us/academic"

# Columns to extract from dbNSFP, mapped to output column names
COLUMN_RENAME_MAP = {
    "Orphanet_disorder_id": "DBNSFP_orphanet_disorder_id",
    "Orphanet_disorder": "DBNSFP_orphanet_disorder",
    "Orphanet_association_type": "DBNSFP_orphanet_association_type",
    "HPO_id": "DBNSFP_hpo_id",
    "HPO_name": "DBNSFP_hpo_name",
    "Function_description": "DBNSFP_function_description",
    "Disease_description": "DBNSFP_disease_description",
    "Pathway(KEGG)_full": "DBNSFP_pathway_kegg",
    "Pathway(Uniprot)": "DBNSFP_pathway_uniprot",
    "P(rec)": "DBNSFP_p_rec",
    "Known_rec_info": "DBNSFP_known_rec_info",
    "Essential_gene": "DBNSFP_essential_gene",
    "MGI_mouse_phenotype": "DBNSFP_mgi_mouse_phenotype",
    "ZFIN_zebrafish_phenotype_tag": "DBNSFP_zfin_zebrafish_phenotype",
}

# Columns used to determine if a gene is disease-associated
DISEASE_FILTER_COLUMNS = ["Orphanet_disorder", "Disease_description", "HPO_id"]


def is_non_empty(value):
    """Check if a value is non-empty (not NaN, not '.', not empty string)."""
    if pd.isna(value):
        return False
    return str(value).strip() not in ("", ".")


@cache_data_table
def get_dbnsfp_gene_table():
    dbnsfp_key = os.environ.get("DBNSFP_KEY", "")
    if not dbnsfp_key:
        raise ValueError("DBNSFP_KEY environment variable is not set")

    url = f"{DBNSFP_BASE_URL}/{dbnsfp_key}/dbNSFP5.3_gene.gz"
    print(f"Downloading {url}")
    response = requests.get(url)
    response.raise_for_status()
    df = pd.read_table(io.BytesIO(response.content), compression="gzip", dtype=str)
    print(f"Read {len(df):,d} rows with {len(df.columns)} columns")

    # Filter to disease-associated genes
    disease_mask = df[DISEASE_FILTER_COLUMNS].apply(lambda col: col.apply(is_non_empty)).any(axis=1)
    df = df[disease_mask]
    print(f"Kept {len(df):,d} disease-associated genes")

    # Handle multi-ENSG rows: explode on ";" and deduplicate
    df["Ensembl_gene"] = df["Ensembl_gene"].str.split(";")
    df = df.explode("Ensembl_gene")
    df["Ensembl_gene"] = df["Ensembl_gene"].str.strip()
    df = df[df["Ensembl_gene"].notna() & (df["Ensembl_gene"] != "") & (df["Ensembl_gene"] != ".")]
    df = df.drop_duplicates(subset=["Ensembl_gene"], keep="first")

    # Select and rename columns
    columns_to_keep = ["Ensembl_gene"] + list(COLUMN_RENAME_MAP.keys())
    df = df[[c for c in columns_to_keep if c in df.columns]]
    df = df.rename(columns=COLUMN_RENAME_MAP)

    # Replace "." null values with empty string
    for col in df.columns:
        if col == "Ensembl_gene":
            continue
        df[col] = df[col].apply(lambda x: "" if pd.isna(x) or str(x).strip() == "." else str(x).strip())

    # Normalize semicolon separators to "; " (with space) for consistency
    string_cols = [c for c in df.columns if c != "Ensembl_gene" and c != "DBNSFP_p_rec"]
    for col in string_cols:
        df[col] = df[col].str.replace(r";\s*", "; ", regex=True)

    df = df.rename(columns={"Ensembl_gene": "gene_id"})

    print(f"Returning {len(df):,d} rows")
    return df


if __name__ == "__main__":
    df = get_dbnsfp_gene_table()
    print(df.shape)
    print(df.head())
