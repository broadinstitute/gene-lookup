"""Download the gene dosage sensitivity scores from Collins et al. 2022 (https://doi.org/10.1016/j.cell.2022.06.036).

pHaplo and pTriplo are the predicted probabilities that losing one copy (pHaplo) or gaining a copy (pTriplo)
of a gene causes disease, estimated from rare copy number variants across 54 disorders. The paper calls
genes with pHaplo >= 0.86 haploinsufficient and genes with pTriplo >= 0.94 triplosensitive (these cutoffs
are in CONSTRAINT_THRESHOLDS in website/global_constants.py, which the website uses for highlighting).

The file covers 18,641 autosomal protein-coding genes, identified only by their GENCODE v19 (GRCh37) gene
symbols, so they have to be mapped to current Ensembl gene ids through HGNC's previous symbols and aliases.
It was published once on Zenodo (record 6347673, CC BY 4.0) and is not updated.

Example rows:

    #gene     pHaplo            pTriplo
    CACNA1C   0.99898184581082  1
    ZNF462    1                 0.987995995573708
"""

import io

import pandas as pd
import requests

from annotation_utils.cache_utils import cache_data_table, record_source_version, return_stale_cached_table_or_raise

COLLINS_DOSAGE_SENSITIVITY_SCORES_URL = (
    "https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1")


def parse_collins_dosage_sensitivity_scores(file_contents):
    """Parse the gzipped Collins et al. 2022 scores table.

    Args:
        file_contents (bytes): the downloaded .tsv.gz file

    Returns:
        pandas.DataFrame: columns gene_symbol, pHaplo, pTriplo
    """
    df = pd.read_table(io.BytesIO(file_contents), compression="gzip")
    df = df.rename(columns={"#gene": "gene_symbol"})
    if df["gene_symbol"].duplicated().any():
        raise ValueError(f"Duplicate gene symbols in the Collins et al. 2022 scores: {df[df['gene_symbol'].duplicated()]['gene_symbol'].tolist()}")
    return df[["gene_symbol", "pHaplo", "pTriplo"]]


@cache_data_table
def get_collins_dosage_sensitivity_scores():
    """Download the Collins et al. 2022 pHaplo and pTriplo scores.

    Returns:
        pandas.DataFrame: one row per gene symbol, with columns gene_symbol, pHaplo, pTriplo
    """
    try:
        response = requests.get(COLLINS_DOSAGE_SENSITIVITY_SCORES_URL, timeout=120)
        response.raise_for_status()
        df = parse_collins_dosage_sensitivity_scores(response.content)
    except Exception as e:
        return return_stale_cached_table_or_raise("get_collins_dosage_sensitivity_scores", "Collins et al. 2022", e)

    # The scores were published once and never revised, so the publication date is the version.
    record_source_version("collins_2022", date_label="published: 2022-03-11")

    return df


if __name__ == "__main__":
    print(get_collins_dosage_sensitivity_scores())
