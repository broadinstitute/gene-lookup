from annotation_utils.cache_utils import cache_data_table
import collections
import requests
from pprint import pprint
import pandas as pd

"""
{'acmg_secondary_findings': [],
'chr': '1',
'clingen_diseases': [],
'clingen_gene_ds': {'entrez_id': 375790,
                    'hgnc_symbol': 'AGRN',
                    'hi_disease_name': None,
                    'hi_morbid_id': None,
                    'hi_pubmed_ids': [],
                    'hi_value': 30,
                    'review_date': '2016-08-22',
                    'ts_disease_name': None,
                    'ts_morbid_id': None,
                    'ts_pubmed_ids': [],
                    'ts_value': 0},
'current_hgnc_id': 329,
'current_hgnc_symbol': 'AGRN',
'end': 1056119,
'ensembl_gene_ensg': 'ENSG00000188157',
'ensembl_hgnc_id': 329,
'ensembl_hgnc_symbol': 'AGRN',
'entrez_id': 375790,
'g2p': None,
'g2p_types': [],
'genccs': [{'classification_name': 'Strong',
            'disease_name': 'congenital myasthenic '
                            'syndrome 8',
            'ensembl_gene_ensg': 'ENSG00000188157',
            'ensembl_hgnc_symbol': 'AGRN',
            'evaluated': '2020-10-26',
            'gencc_acc': 'GENCC_000111-HGNC_329-OMIM_615120-HP_0000007-GENCC_100002',
            'inheritance_mode': 'Autosomal recessive '
                                'inheritance',
            'mondo_id': 14052,
            'submitter_name': 'PanelApp Australia'},
            {'classification_name': 'Strong',
            'disease_name': 'congenital myasthenic '
                            'syndrome 8',
            'ensembl_gene_ensg': 'ENSG00000188157',
            'ensembl_hgnc_symbol': 'AGRN',
            'evaluated': '2021-01-29',
            'gencc_acc': 'GENCC_000104-HGNC_329-OMIM_615120-HP_0000007-GENCC_100002',
            'inheritance_mode': 'Autosomal recessive '
                                'inheritance',
            'mondo_id': 14052,
            'submitter_name': 'Genomics England '
                                'PanelApp'},
            {'classification_name': 'Supportive',
            'disease_name': 'postsynaptic congenital '
                            'myasthenic syndrome',
            'ensembl_gene_ensg': 'ENSG00000188157',
            'ensembl_hgnc_symbol': 'AGRN',
            'evaluated': '2021-09-14',
            'gencc_acc': 'GENCC_000110-HGNC_329-Orphanet_98913-HP_0000007-GENCC_100009',
            'inheritance_mode': 'Autosomal recessive '
                                'inheritance',
            'mondo_id': 20344,
            'submitter_name': 'Orphanet'},
            {'classification_name': 'Supportive',
            'disease_name': 'obsolete presynaptic '
                            'congenital myasthenic '
                            'syndrome',
            'ensembl_gene_ensg': 'ENSG00000188157',
            'ensembl_hgnc_symbol': 'AGRN',
            'evaluated': '2021-09-14',
            'gencc_acc': 'GENCC_000110-HGNC_329-Orphanet_98914-HP_0000006-GENCC_100009',
            'inheritance_mode': 'Autosomal dominant '
                                'inheritance',
            'mondo_id': 20345,
            'submitter_name': 'Orphanet'},
            {'classification_name': 'Strong',
            'disease_name': 'congenital myasthenic '
                            'syndrome 8',
            'ensembl_gene_ensg': 'ENSG00000188157',
            'ensembl_hgnc_symbol': 'AGRN',
            'evaluated': '2018-12-15',
            'gencc_acc': 'GENCC_000106-HGNC_329-OMIM_615120-HP_0000007-GENCC_100002',
            'inheritance_mode': 'Autosomal recessive '
                                'inheritance',
            'mondo_id': 14052,
            'submitter_name': 'Labcorp Genetics '
                                '(formerly Invitae)'}],
'hgnc_description': 'agrin',
'hi_index': 71.12874595384775,
'is_protein_coding': True,
'loeuf': 0.73,
'omim_ids': [103320],
'omim_morbid_diseases': [{'association_type': 3,
                            'disease_name': 'Myasthenic '
                                            'syndrome, '
                                            'congenital, '
                                            '8, with '
                                            'pre- and '
                                            'postsynaptic '
                                            'defects',
                            'inheritance_mode': 'Autosomal '
                                                'recessive',
                            'omim_morbid_id': 615120}],
'p_hi': '0.784779631672505',
'p_li': 1.5815e-20,
'p_ts': '0.296900840763816',
'public_snv_count': 5,
's_het': '0.014729',
'start': 1020120,
'ucsc_acc': 'uc001ack.3',
'uniprot_accs': ['O00468']},
"""


def rename_inheritance_mode(inheritance_mode):
    if inheritance_mode == "Autosomal dominant inheritance":
        return "AD"
    elif inheritance_mode == "Autosomal recessive inheritance":
        return "AR"
    elif inheritance_mode == "X-linked inheritance" or inheritance_mode == "X-linked recessive inheritance":
        return "XR"
    elif inheritance_mode == "Semidominant inheritance":
        return "SD"
    elif inheritance_mode == "Mitochondrial inheritance":
        return "MITO"
    else:
        return ""


@cache_data_table
def get_decipher_gene_table():
    url = "https://www.deciphergenomics.org/data/genes"
    response = requests.get(url, headers={
        "accept": "application/json",
        "x-requested-with": "XMLHttpRequest",
    })

    if not response.ok:
        raise ValueError(f"Failed to get data from {url}")

    data = response.json()

    genes = data.get("content").get("genes")

    if not genes:
        raise ValueError(f"Failed to get genes from {url}: {data}")

    #inheritance_mode_counter = collections.Counter()
    output_records = []
    for gene in genes:
        if "genccs" not in gene:
            continue

        gene_id_to_disease_name = collections.defaultdict(set)

        for gencc in gene.get("genccs", []):
            gene_id = gencc.get("ensembl_gene_ensg")
            disease_name = gencc.get("disease_name")
            inheritance_mode = rename_inheritance_mode(gencc.get("inheritance_mode", ""))

            if gene_id and disease_name:
                gene_id_to_disease_name[gene_id].add((inheritance_mode, disease_name))
                #if inheritance_mode:
                #    inheritance_mode_counter[inheritance_mode] += 1

        for gene_id, disease_name_and_inheritance_mode in gene_id_to_disease_name.items():
            sorted_disease_name_and_inheritance_mode = sorted(disease_name_and_inheritance_mode)
            inheritance_modes = "; ".join([x[0] for x in sorted_disease_name_and_inheritance_mode])
            disease_names = "; ".join([x[1] for x in sorted_disease_name_and_inheritance_mode])

            output_records.append({
                "gene_id": gene_id,
                "inheritance_modes": inheritance_modes,
                "disease_names": disease_names,
            })

        
    return pd.DataFrame(output_records)
        

if __name__ == "__main__":
    df = get_decipher_gene_table()

    print(df)


