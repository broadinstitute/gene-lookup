import hail as hl
import hailtop.fs as hfs
import os
import pandas as pd
from annotation_utils.cache_utils import cache_data_table

"""
----------------------------------------
Global fields:
    'clinvar_release_date': str
    'mane_select_version': str
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'clinvar_variation_id': str
    'rsid': str
    'review_status': str
    'gold_stars': int32
    'clinical_significance': str
    'last_evaluated': str
    'submissions': array<struct {
        id: str,
        submitter_name: str,
        clinical_significance: str,
        last_evaluated: str,
        review_status: str,
        conditions: array<struct {
            name: str,
            medgen_id: str
        }>
    }>
    'variant_id': str
    'reference_genome': str
    'chrom': str
    'pos': int32
    'ref': str
    'alt': str
    'transcript_consequences': array<struct {
        biotype: str,
        consequence_terms: array<str>,
        domains: set<str>,
        gene_id: str,
        gene_symbol: str,
        hgvsc: str,
        hgvsp: str,
        is_canonical: bool,
        lof_filter: str,
        lof_flags: str,
        lof: str,
        major_consequence: str,
        transcript_id: str,
        polyphen_prediction: str,
        sift_prediction: str,
        transcript_version: str,
        gene_version: str,
        is_mane_select: bool,
        is_mane_select_version: bool,
        refseq_id: str,
        refseq_version: str
    }>
    'gnomad': struct {
        exome: struct {
            filters: set<str>,
            ac: int32,
            an: int32
        },
        genome: struct {
            filters: set<str>,
            ac: int32,
            an: int32
        }
    }
    'in_gnomad': bool
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------


collections.Counter(ht.clinical_significance.collect())

Counter({'Uncertain significance': 1778465,
         'Likely benign': 986724,
         'Benign': 209170,
         'Pathogenic': 170195,
         'Conflicting classifications of pathogenicity': 144798,
         'Likely pathogenic': 100497,
         'Benign/Likely benign': 53343,
         'Pathogenic/Likely pathogenic': 34349,
         'not provided': 10875,
         'drug response': 1849,
         'other': 1547,
         'risk factor': 381,
         'association': 339,
         'Uncertain significance/Uncertain risk allele': 148,
         'no classifications from unflagged records': 144,
         'Affects': 135,
         'Likely risk allele': 94,
         'Pathogenic; other': 73,
         'Pathogenic; drug response': 55,
         'protective': 38,
         'Uncertain risk allele': 33,
         ...

"""

@cache_data_table
def get_clinvar_gene_disease_table():
    # ht.clinvar_release_date.collect()

    # Filter to transcripts that are MANE Select or MANE Plus Clinical
    # List downloaded from http://tark.ensembl.org/web/manelist/
    df_mane = pd.read_csv(os.path.join(os.path.dirname(__file__), "data", "MANE_select_and_plus_ENSG_ids.csv"))
    mane_transcripts = [transcript_id.split(".")[0] for transcript_id in df_mane["ENST_ids"]]


    ht = hl.read_table("gs://gnomad-v4-data-pipeline/output/clinvar/clinvar_grch38_annotated_2.ht")

    # check if clinical_significance string (when converted to lower case) contains "pathogenic" but not "pathogenicity"
    ht = ht.filter(hl.str(ht.clinical_significance).lower().contains("pathogenic") & ~hl.str(ht.clinical_significance).lower().contains("pathogenicity"), keep=True)

    # get the set of ENSG gene ids from the transcript_consequences field using 
    ht = ht.annotate(gene_id = hl.array(hl.set(ht.transcript_consequences.filter(lambda x: x.gene_id.startswith("ENSG")).map(lambda x: x.gene_id))))  #  & hl.set(mane_transcripts).contains(x.transcript_id)

    # filter out rows where gene_ids is empty
    ht = ht.filter(hl.len(ht.gene_id) > 0, keep=True)
    ht = ht.annotate(phenotypes=hl.array(hl.sorted(hl.set(ht.submissions.map(lambda x: hl.str(", ").join(x.conditions.map(lambda y: y.name)))))))
    # filter out rows where phenotypes is empty
    ht = ht.filter(hl.len(ht.phenotypes) > 0, keep=True)
    ht = ht.annotate(phenotypes = hl.str(", ").join(ht.phenotypes))
    ht = ht.filter(ht.phenotypes != "not provided", keep=True)

    ht = ht.annotate(major_consequences = hl.str(", ").join(hl.array(hl.sorted(hl.set(ht.transcript_consequences.map(lambda x: x.major_consequence))))))

    ht = ht.key_by()

    ht = ht.select(
        "gene_id",
        "phenotypes",
        "major_consequences",
        "clinical_significance",
        #"last_evaluated",
        #"review_status",
        "gold_stars",
    )
    ht = ht.explode("gene_id")
    df = ht.to_pandas()

    df.drop_duplicates(inplace=True)
    
    # group by gene_id and combine the other fields using ; as a separator
    df = df.groupby("gene_id").agg({
        "phenotypes": lambda x: ", ".join(list(sorted(set(", ".join(x).split(", "))))[:5]),
        "clinical_significance": lambda x: ", ".join(list(sorted(set(", ".join(x).split(", "))))[:5]),
        "gold_stars": lambda x: max(x),
        "major_consequences": lambda x: ", ".join(list(sorted(set(", ".join(x).split(", "))))[:5]),
    })

    df.reset_index(inplace=True)

    return df


def simplify_clinical_significance_expr(clinsig_expr):
    return (
        hl.case()
          .when(clinsig_expr.contains("Pathogenic"), "Pathogenic")
          .when(clinsig_expr.contains("Likely pathogenic"), "Likely pathogenic")
          .when(clinsig_expr.contains("Benign"), "Benign")
          .when(clinsig_expr.contains("Likely benign"), "Likely benign")
          .when(clinsig_expr.contains("Uncertain significance"), "Uncertain significance")
          .default(clinsig_expr)
    )

def export_clinvar_vcf(only_pathogenic=True, include_phenotypes=True):
    """Generate a VCF from the latest clinvar hail table in the gnomAD bucket"""
    hl.init(driver_memory="highmem", idempotent=True)
    ht = hl.read_table("gs://gnomad-v4-data-pipeline/output/clinvar/clinvar_grch38_annotated_2.ht")
    ht = ht.key_by(ht.locus, ht.alleles)

    clinvar_release_date = ht.clinvar_release_date.collect()
    clinvar_release_date = clinvar_release_date[0].replace("-", "_")

    print(f"Exporting ClinVar VCF for release date: {clinvar_release_date}")
    # check if clinical_significance string (when converted to lower case) contains "pathogenic" but not "pathogenicity"
    # this excludes 'Conflicting classifications of pathogenicity' as well as 'Uncertain significance'
    # but keeps 'Pathogenic', 'Likely pathogenic', and the many other combinations of these with other labels
    if only_pathogenic:
        ht = ht.filter(hl.str(ht.clinical_significance).lower().contains("pathogenic") & ~hl.str(ht.clinical_significance).lower().contains("pathogenicity"), keep=True)
        ht = ht.annotate(clinical_significance = simplify_clinical_significance_expr(ht.clinical_significance))

    #ht = ht.annotate(major_consequences = hl.str(",").join(hl.array(hl.sorted(hl.set(ht.transcript_consequences.map(lambda x: x.major_consequence))))))
    if include_phenotypes:
        ht = ht.annotate(phenotypes = ht.submissions
            .flatmap(lambda x: x.conditions.map(lambda y: y.name))
            .filter(lambda p: (p.lower() != "not provided") & (p.lower() != "not specified"))
        )
        ht = ht.annotate(phenotypes = hl.str(", ").join(hl.sorted(hl.set(ht.phenotypes))).replace(" ", "_"))

    
    ht = ht.transmute(info=hl.struct(
        clinsig=hl.str(ht.clinical_significance).replace(" ", "_"),
        stars=ht.gold_stars,
        clinvarid=ht.clinvar_variation_id,
        #consequences=ht.major_consequences,
        #review_status=ht.review_status,  # this is equivalent to stars:  practice guideline (4 stars), reviewed by expert panel (3 stars), criteria provided, multiple submitters (2 stars), criteria provided, single submitter (1 star), no assertion criteria provided (0 stars)
        #last_evaluated=ht.last_evaluated,
        #in_gnomad=ht.in_gnomad,
    ))

    if include_phenotypes:
        ht = ht.transmute(info=ht.info.annotate(phenotypes=ht.phenotypes))

    ht = ht.select('info')
    ht = ht.checkpoint(f"gs://bw2-delete-after-5-days/clinvar_{clinvar_release_date}_checkpoint.ht", overwrite=True)
    local_filename = f"clinvar_{clinvar_release_date}.vcf.bgz"
    temp_bucket_path = f"gs://bw2-delete-after-5-days/{local_filename}"

    hl.export_vcf(ht, temp_bucket_path)
    hfs.copy(temp_bucket_path, f"./{local_filename}")
    hfs.remove(temp_bucket_path)

    return local_filename


if __name__ == "__main__":
    #df = get_clinvar_gene_disease_table()
    #print(df)
    local_filename = export_clinvar_vcf(only_pathogenic=True)
    print(f"Exported ClinVar VCF to {local_filename}")