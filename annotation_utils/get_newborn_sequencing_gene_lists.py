"""Which genomic newborn screening programs include each gene, and the gene's RUSP category.

Source: Minten et al. 2025, "Data-driven consideration of genetic disorders for global genomic newborn screening
programs", Genetics in Medicine 27:101443. The authors collected the gene lists of 27 newborn sequencing programs
(eg. BabySeq, GUARDIAN, the Generation Study, BabyScreen+), linked each gene to one disorder, and recorded
whether that disorder is a core or secondary condition of the US Recommended Uniform Screening Panel (RUSP).
Their data is at https://github.com/tjminten/gim_genelist (Data_Repository.xlsx).

The data is a one-time snapshot, so the table derived from it is stored in the repository
(NEWBORN_SEQUENCING_GENE_LISTS_PATH) rather than downloaded on every run. To regenerate it, run this module
(needs the openpyxl package to read the Excel file).
"""

import io
import re

import pandas as pd
import requests

NEWBORN_SEQUENCING_GENE_LISTS_PATH = "annotation_utils/data/newborn_sequencing_gene_lists_from_Minten_2025.tsv"
MINTEN_2025_DATA_REPOSITORY_URL = (
    "https://raw.githubusercontent.com/tjminten/gim_genelist/e6a05dafa1c4bb4de27c40f41952c37e84125fa8/Data_Repository.xlsx")


def convert_minten_data_repository_to_gene_table(df_repository, df_codebook):
    """Summarize the Minten et al. 2025 data repository as one row per gene.

    Args:
        df_repository (pandas.DataFrame): the "Repository" sheet: one row per gene, with a 1 in each scr_<program>
            column whose program includes the gene, and the RUSP category in the rusp column
        df_codebook (pandas.DataFrame): the "Codebook" sheet, read without a header: column 0 is the variable name
            and column 1 its description, eg. "Inclusion on BabySeq gene list (1 if gene is present)"

    Returns:
        pandas.DataFrame: columns gene_symbol, rusp_category ("Core", "Secondary" or ""), program_count, and
            programs (the names of the programs that include the gene, joined with ", ")
    """
    program_columns = [column for column in df_repository.columns if column.startswith("scr_")]
    variable_descriptions = dict(zip(df_codebook[0].astype(str), df_codebook[1].astype(str)))
    program_names = {}
    for column in program_columns:
        match = re.match(r"Inclusion on (.+?) gene list", variable_descriptions.get(column, ""))
        if not match:
            raise ValueError(f"No program name for {column} in the codebook")
        program_names[column] = match.group(1)

    included = df_repository[program_columns].fillna(0).astype(int) == 1
    program_counts = included.sum(axis=1)
    if not (program_counts == df_repository["programs_screened"]).all():
        raise ValueError("programs_screened doesn't match the number of programs marked in the scr_ columns")

    return pd.DataFrame({
        # The repository writes open reading frame symbols in upper case (C19ORF12); HGNC writes C19orf12
        "gene_symbol": df_repository["gene"].str.replace(r"^(C(?:\d+|X|Y))ORF(\d+)$", r"\1orf\2", regex=True),
        # ANALYSIS_OK[imputation]: "" means the gene's disorder is not a RUSP condition
        "rusp_category": df_repository["rusp"].fillna(""),
        "program_count": program_counts,
        "programs": included.apply(lambda row: ", ".join(program_names[column] for column in program_columns if row[column]), axis=1),
    })


def read_newborn_sequencing_gene_lists():
    """Return the table stored at NEWBORN_SEQUENCING_GENE_LISTS_PATH (see convert_minten_data_repository_to_gene_table)."""
    return pd.read_table(NEWBORN_SEQUENCING_GENE_LISTS_PATH, keep_default_na=False)


if __name__ == "__main__":
    response = requests.get(MINTEN_2025_DATA_REPOSITORY_URL, timeout=120)
    response.raise_for_status()
    workbook = pd.ExcelFile(io.BytesIO(response.content))
    df_gene_lists = convert_minten_data_repository_to_gene_table(
        workbook.parse("Repository"), workbook.parse("Codebook", header=None))
    df_gene_lists.to_csv(NEWBORN_SEQUENCING_GENE_LISTS_PATH, sep="\t", index=False)
    print(f"Wrote {len(df_gene_lists):,d} genes to {NEWBORN_SEQUENCING_GENE_LISTS_PATH}")
