import pandas as pd
from annotation_utils.cache_utils import cache_data_table
from annotation_utils.get_hgnc_table import get_hgnc_table

"""Retrieve GenCC submissions

uuid                                   GENCC_000101-HGNC_10896-OMIM_182212-HP_0000006...
** gene_curie                                                                    HGNC:10896
gene_symbol                                                                          SKI
** disease_curie                                                              MONDO:0008426
** disease_title                                               Shprintzen-Goldberg syndrome
disease_original_curie                                                       OMIM:182212
disease_original_title                                      Shprintzen-Goldberg syndrome
classification_curie                                                        GENCC:100001
** classification_title                                                          Definitive
moi_curie                                                                     HP:0000006
** moi_title                                                             Autosomal dominant
submitter_curie                                                             GENCC:000101
submitter_title                                                           Ambry Genetics
submitted_as_hgnc_id                                                          HGNC:10896
submitted_as_hgnc_symbol                                                             SKI
submitted_as_disease_id                                                      OMIM:182212
submitted_as_disease_name                                   Shprintzen-Goldberg syndrome
submitted_as_moi_id                                                           HP:0000006
submitted_as_moi_name                                     Autosomal dominant inheritance
submitted_as_submitter_id                                                   GENCC:000101
submitted_as_submitter_name                                               Ambry Genetics
submitted_as_classification_id                                              GENCC:100001
submitted_as_classification_name                                              Definitive
submitted_as_date                                                    2018-03-30 13:31:56
submitted_as_public_report_url                                                       NaN
submitted_as_notes                                                                   NaN
submitted_as_pmids                                                                   NaN
submitted_as_assertion_criteria_url                                       PMID: 28106320
submitted_as_submission_id                                                          1034
submitted_run_date                                                            2020-12-24
"""

def normalize_nulls(x):
    if pd.isna(x):
        return ""
    return x



@cache_data_table
def get_gencc_table():
    df_hgnc = get_hgnc_table()
    HGNC_to_ENSG_map = dict(zip(df_hgnc["HGNC ID"], df_hgnc["Ensembl gene ID"]))

    df = pd.read_table( "https://search.thegencc.org/download/action/submissions-export-tsv")

    df["gene_id"] = df["gene_curie"].map(HGNC_to_ENSG_map)
    df["hgnc_gene_id"] = df["gene_curie"]
    df.rename(columns={
        "disease_curie": "disease_id",
        "disease_title": "disease_name",
        "classification_title": "classification",
        "moi_title": "inheritance",
    }, inplace=True)

    df = df[[
        "gene_id",
        "hgnc_gene_id",
        "disease_id",
        "disease_name",
        "classification",
        "inheritance",
    ]]

    """
    Autosomal recessive     10132
    Autosomal dominant       8913
    X-linked                 1004
    Unknown                   928
    Semidominant              139
    Mitochondrial              98
    Y-linked inheritance        2
    X-linked recessive          2
    """
    # rename using dictionary
    df["inheritance"] = df["inheritance"].map({
        "Autosomal recessive": "AR",
        "Autosomal dominant": "AD",
        "Semidominant": "SD",
        "Mitochondrial": "MITO",
        "X-linked recessive": "XR",
        "X-linked": "XR",
        "Y-linked inheritance": "YL",
    })

    return df


if __name__ == "__main__":
    df = get_gencc_table()
    print(df)
    print(df["inheritance"].value_counts())
