import pandas as pd
from annotation_utils.cache_utils import cache_data_table
from annotation_utils.get_ensembl_db_info import get_transcript_id_to_gene_id


@cache_data_table
def _get_gnomAD_v2_constraint():
    transcript_id_to_gene_id = get_transcript_id_to_gene_id(database="homo_sapiens_core_74_37")

    # gnomAD v2.1.1 lof metrics
    df = pd.read_table("https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression="gzip")
    df["gene_id"] = df["transcript"].map(transcript_id_to_gene_id)
    before = len(df)
    examples = df[df["gene_id"].isna()]["transcript"][:10]
    df = df[df["gene_id"].notna()]
    if before - len(df) > 0:
        print(f"WARNING: Dropped {before - len(df):,d} out of {before:,d} ({((before - len(df)) / before) * 100:.2f}%) rows from gnomAD v2.1.1 lof metrics because the transcript id was not found in the Ensembl database. Examples: {', '.join(examples)}")

    df = df[df["gene_id"].str.startswith("ENSG")]

    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    print(df.iloc[0])

    df = df[[
        "gene_id",
        "pLI",
    ]]

    df = df.rename(columns={
        "pLI": "pLI_v2",
    })

    df = df[df["pLI_v2"].notna()]
    df = df.sort_values(by="pLI_v2", ascending=False)
    df = df.drop_duplicates(subset=["gene_id"], keep="first")

    """
    Columns:
        gene                            MED13
        transcript            ENST00000397786
        obs_mis                           871
        exp_mis                        1117.8
        oe_mis                        0.77921
        mu_mis                       0.000056
        possible_mis                    14195
        obs_mis_pphen                   314.0
        exp_mis_pphen                  529.75
        oe_mis_pphen                  0.59273
        possible_mis_pphen             6708.0
        obs_syn                           422
        exp_syn                        387.53
        oe_syn                          1.089
        mu_syn                       0.000019
        possible_syn                     4248
        obs_lof                           0.0
        mu_lof                       0.000005
        possible_lof                   1257.0
        exp_lof                        98.429
        pLI                               1.0
        pNull                             0.0
        pRec                              0.0
        oe_lof                            0.0
        oe_syn_lower                    1.005
        oe_syn_upper                     1.18
        oe_mis_lower                    0.736
        oe_mis_upper                    0.824
        oe_lof_lower                      0.0
        oe_lof_upper                     0.03
        constraint_flag                   NaN
        syn_z                         -1.3765
        mis_z                          2.6232
        lof_z                          9.1935
        oe_lof_upper_rank                 0.0
        oe_lof_upper_bin                  0.0
        oe_lof_upper_bin_6                0.0
        n_sites                           2.0
        classic_caf                  0.000012
        max_af                       0.000008
        no_lofs                      124782.0
        obs_het_lof                       3.0
        obs_hom_lof                       0.0
        defined                      124785.0
        p                            0.000012
        exp_hom_lof                  0.000018
        classic_caf_afr                   0.0
        classic_caf_amr                   0.0
        classic_caf_asj                   0.0
        classic_caf_eas                   0.0
        classic_caf_fin              0.000093
        classic_caf_nfe              0.000009
        classic_caf_oth                   0.0
        classic_caf_sas                   0.0
        p_afr                             0.0
        p_amr                             0.0
        p_asj                             0.0
        p_eas                             0.0
        p_fin                        0.000093
        p_nfe                        0.000009
        p_oth                             0.0
        p_sas                             0.0
        transcript_type        protein_coding
        gene_id               ENSG00000108510
        transcript_level                    2
        cds_length                       6522
        num_coding_exons                   30
        gene_type              protein_coding
        gene_length                    122678
        exac_pLI                          1.0
        exac_obs_lof                      0.0
        exac_exp_lof                   64.393
        exac_oe_lof                       0.0
        brain_expression                  NaN
        chromosome                         17
        start_position               60019966
        end_position                 60142643    
    """

    return df


@cache_data_table
def _get_gnomAD_v4_constraint():
     # gnomAD v4.1 constraint metrics
    df = pd.read_table("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv")
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    print(df2.iloc[0])

    df = df[df["gene_id"].str.startswith("ENSG")]

    df = df[[
        "gene_id",
        "lof.pLI",
        "lof.oe",
        "lof.oe_ci.lower",
        "lof.oe_ci.upper",
        "mis.oe",
        "mis.oe_ci.lower",
        "mis.oe_ci.upper",
    ]]

    df = df.rename(columns={
        "lof.pLI": "pLI_v4",
        "lof.oe": "lof_oe_v4",
        "lof.oe_ci.lower": "lof_oe_ci_lower_v4",
        "lof.oe_ci.upper": "lof_oe_ci_upper_v4",
        "mis.oe": "mis_oe_v4",
        "mis.oe_ci.lower": "mis_oe_ci_lower_v4",
        "mis.oe_ci.upper": "mis_oe_ci_upper_v4",
    })


    df = df[df["pLI_v4"].notna()]
    df = df.sort_values(by=["pLI_v4", "lof_oe_v4"], ascending=[False, True])
    df = df.drop_duplicates(subset=["gene_id"], keep="first")

    """
    Columns:
        gene                                 A1BG
        gene_id                                 1
        transcript                    NM_130786.4
        canonical                            True
        mane_select                          True
        lof_hc_lc.obs                        45.0
        lof_hc_lc.exp                      43.048
        lof_hc_lc.possible                  193.0
        lof_hc_lc.oe                       1.0454
        lof_hc_lc.mu                     0.000001
        lof_hc_lc.pLI                         0.0
        lof_hc_lc.pNull                   0.84915
        lof_hc_lc.pRec                    0.15085
        lof.obs                              45.0
        lof.exp                            43.048
        lof.possible                        193.0
        lof.oe                             1.0454
        lof.mu                           0.000001
        lof.pLI                               0.0
        lof.pNull                         0.84295
        lof.pRec                          0.15705
        lof.oe_ci.lower                     0.822
        lof.oe_ci.upper                      1.34
        lof.oe_ci.upper_rank                  NaN
        lof.oe_ci.upper_bin_decile            NaN
        lof.z_raw                        -0.29756
        lof.z_score                      -0.25212
        mis.obs                             707.0
        mis.exp                            647.03
        mis.possible                       2870.0
        mis.oe                             1.0927
        mis.mu                           0.000008
        mis.oe_ci.lower                     1.026
        mis.oe_ci.upper                     1.163
        mis.z_raw                         -2.3574
        mis.z_score                      -0.86092
        mis_pphen.obs                       220.0
        mis_pphen.exp                       190.6
        mis_pphen.possible                  890.0
        mis_pphen.oe                       1.1543
        syn.obs                             316.0
        syn.exp                            295.94
        syn.possible                        994.0
        syn.oe                             1.0678
        syn.mu                           0.000003
        syn.oe_ci.lower                     0.973
        syn.oe_ci.upper                     1.172
        syn.z_raw                          -1.166
        syn.z_score                      -0.63549
        constraint_flags                       []
        level                                 NaN
        transcript_type                       NaN
        chromosome                            NaN
        cds_length                            NaN
        num_coding_exons                      NaN
    """

    return df


@cache_data_table
def _get_s_het_scores():
    """Download s_het estimates from Zeng et al. 2024 (GeneBayes).

    Source: https://zenodo.org/records/7939768
    Paper: https://doi.org/10.1038/s41588-024-01820-9
    Scores >0.1 indicate high likelihood of extreme selection.
    """
    df = pd.read_table("https://zenodo.org/api/records/7939768/files/s_het_estimates.genebayes.tsv/content")
    df = df[df["ensg"].str.startswith("ENSG")]
    df = df[["ensg", "post_mean"]].rename(columns={"ensg": "gene_id", "post_mean": "s_het"})
    df = df[df["s_het"].notna()]
    df = df.drop_duplicates(subset=["gene_id"], keep="first")
    return df


def get_constraint_scores():
    "Download gnomAD gene constraint tables and s_het scores"

    df = _get_gnomAD_v2_constraint()
    df2 = _get_gnomAD_v4_constraint()
    df3 = _get_s_het_scores()

    df = df.set_index("gene_id").join(df2.set_index("gene_id"), how="outer")
    df = df.join(df3.set_index("gene_id"), how="outer").reset_index()

    return df


if __name__ == "__main__":
    df = get_constraint_scores()
    print(df.head())