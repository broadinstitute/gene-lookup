from annotation_utils.cache_utils import cache_data_table
import pandas as pd

MANE_SUMMARY_TABLE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz"


@cache_data_table
def get_MANE_ensembl_transcript_table(mane_summary_table_url=MANE_SUMMARY_TABLE_URL):
    """Download the MANE summary table and return it as a pandas DataFrame."""
    return pd.read_table(mane_summary_table_url)


"""
Columns:

0: #NCBI_GeneID
1: Ensembl_Gene
2: HGNC_ID
3: symbol
4: name
5: RefSeq_nuc
6: RefSeq_prot
7: Ensembl_nuc
8: Ensembl_prot
9: MANE_status
10: GRCh38_chr
11: chr_start
12: chr_end
13: chr_strand

(19404, 14)

Example:

#NCBI_GeneID                  GeneID:1
Ensembl_Gene        ENSG00000121410.13
HGNC_ID                         HGNC:5
symbol                            A1BG
name            alpha-1-B glycoprotein
RefSeq_nuc                 NM_130786.4
RefSeq_prot                NP_570602.2
Ensembl_nuc          ENST00000263100.8
Ensembl_prot         ENSP00000263100.2
MANE_status                MANE Select
GRCh38_chr                NC_000019.10
chr_start                     58345183
chr_end                       58353492
chr_strand                           -
"""