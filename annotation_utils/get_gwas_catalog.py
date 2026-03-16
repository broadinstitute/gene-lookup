# download the GWAS catalog and parse gene disease relationships for the subset of records where the MONDO term is a rare disease term

import io
import pandas as pd
import requests
import os
import zipfile

from annotation_utils.cache_utils import cache_data_table
from annotation_utils.get_mondo_ontology import get_mondo_ontology

GWAS_CATALOG_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"

@cache_data_table
def _download_gwas_catalog():
    """Download the GWAS catalog and return it as a pandas DataFrame."""
    import time as _time

    for attempt in range(1, 4):
        try:
            print(f"Downloading GWAS catalog (attempt {attempt})")
            r = requests.get(GWAS_CATALOG_URL, timeout=120)
            r.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
                return pd.read_table(zf.open(zf.namelist()[0]))
        except Exception as e:
            if attempt < 3:
                print(f"  Failed: {e}. Retrying in 30 seconds...")
                _time.sleep(30)
            else:
                raise


@cache_data_table
def get_gwas_catalog_rare_disease_records():
    """
    Download the GWAS catalog and parse gene disease relationships for the subset of records where the MONDO term is a rare disease term
    """

    """
    #DATE ADDED TO CATALOG                                                2020-08-03
    #PUBMEDID                                                               32352494
    #FIRST AUTHOR                                                              Han X
    #DATE                                                                 2020-04-30
    #JOURNAL                                                         JAMA Ophthalmol
    #LINK                                       www.ncbi.nlm.nih.gov/pubmed/32352494
    #STUDY                         Association of Myopia and Intraocular Pressure...
    #DISEASE/TRAIT                                              Spherical equivalent
    #INITIAL SAMPLE SIZE                        95,827 European ancestry individuals
    #REPLICATION SAMPLE SIZE                                                     NaN
    #REGION                                                                  11p14.1
    CHR_ID                                                                       11
    CHR_POS                                                                30033063
    #REPORTED GENE(S)                                                             NR
    #MAPPED_GENE                                                  KCNA4 - ARL14EP-DT
    UPSTREAM_GENE_ID                                                ENSG00000182255
    DOWNSTREAM_GENE_ID                                              ENSG00000254532
    SNP_GENE_IDS                                                                NaN
    UPSTREAM_GENE_DISTANCE                                                  16033.0
    DOWNSTREAM_GENE_DISTANCE                                                10990.0
    STRONGEST SNP-RISK ALLELE                                           rs7944541-G
    SNPS                                                                  rs7944541
    MERGED                                                                        0
    SNP_ID_CURRENT                                                        7944541.0
    CONTEXT                                                      intergenic_variant
    INTERGENIC                                                                  1.0
    RISK ALLELE FREQUENCY                                                        NR
    P-VALUE                                                                     0.0
    PVALUE_MLOG                                                            15.30103
    P-VALUE (TEXT)                                                              NaN
    OR or BETA                                                                 0.11
    95% CI (TEXT)                                         [0.09-0.13] unit increase
    PLATFORM [SNPS PASSING QC]                                              NR [NR]
    CNV                                                                           N
    MAPPED_TRAIT                                                   refractive error
    MAPPED_TRAIT_URI                   http://purl.obolibrary.org/obo/MONDO_0004892
    STUDY ACCESSION                                                      GCST010378
    GENOTYPING TECHNOLOGY                              Genome-wide genotyping array
    """

    mondo_rare_disease_term_lookup = get_mondo_ontology()

    df_gwas = _download_gwas_catalog()
    df_gwas = df_gwas[df_gwas["MAPPED_TRAIT_URI"].notna()]
    df_gwas["MONDO_ID"] = df_gwas["MAPPED_TRAIT_URI"].apply(os.path.basename).str.replace("_", ":")

    df_gwas = df_gwas[df_gwas["MONDO_ID"].isin(mondo_rare_disease_term_lookup)]
    df_gwas = df_gwas[[
        "MONDO_ID", 
        "CHR_ID",
        "CHR_POS",
        "SNPS",
        "UPSTREAM_GENE_ID",
        "DOWNSTREAM_GENE_ID",
        "UPSTREAM_GENE_DISTANCE",
        "DOWNSTREAM_GENE_DISTANCE",
        "P-VALUE",
        "OR or BETA",
        "95% CI (TEXT)",
    ]]
    df_gwas["MONDO_NAME"] = df_gwas["MONDO_ID"].apply(lambda x: mondo_rare_disease_term_lookup[x]["name"])
    df_gwas["MONDO_CATEGORY"] = df_gwas["MONDO_ID"].apply(lambda x: mondo_rare_disease_term_lookup[x].get("category"))

    # First melt the gene IDs
    df_genes = df_gwas.melt(
        id_vars=["MONDO_ID", "MONDO_NAME", "MONDO_CATEGORY", "CHR_ID", "CHR_POS", "SNPS", "P-VALUE", "OR or BETA", "95% CI (TEXT)", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE"],
        value_vars=["UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID"],
        var_name="GENE_TYPE",
        value_name="GENE_ID"
    )
    df_genes["GENE_TYPE"] = df_genes["GENE_TYPE"].str.replace("_GENE_ID", "")
    df_genes = df_genes[df_genes["GENE_ID"].notna()]

    # Then melt the distances
    df_genes_and_distances = df_genes.melt(
        id_vars=["MONDO_ID", "MONDO_NAME", "MONDO_CATEGORY", "CHR_ID", "CHR_POS", "SNPS", "P-VALUE", "OR or BETA", "95% CI (TEXT)", "GENE_ID", "GENE_TYPE"],
        value_vars=["UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE"],
        var_name="DISTANCE_TYPE",
        value_name="GENE_DISTANCE"
    )

    # Clean up the type columns to match
    df_genes_and_distances["DISTANCE_TYPE"] = df_genes_and_distances["DISTANCE_TYPE"].str.replace("_GENE_DISTANCE", "")
    df_genes_and_distances = df_genes_and_distances[
        df_genes_and_distances["GENE_DISTANCE"].notna() & (df_genes_and_distances["GENE_TYPE"] == df_genes_and_distances["DISTANCE_TYPE"])
    ]

    df_gwas = df_genes_and_distances.drop(columns=["DISTANCE_TYPE"])
    
    return df_gwas
    


