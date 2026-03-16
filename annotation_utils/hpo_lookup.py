import argparse
from annotation_utils.cache_utils import cache_json
import os
import requests
from tqdm import tqdm

HP_OBO_URL = 'http://purl.obolibrary.org/obo/hp.obo'


@cache_json
def download_hpo_obo_file():
    """
    Parse an .obo file which contains a record for each term in the Human Phenotype Ontology

    Args:
        file_iterator: Iterator over lines in the hp.obo file
    Returns:
        dictionary that maps HPO id strings to a record containing
    """

    print("Downloading hp.obo data")
    hpo_id_to_record = {}
    request = requests.get(HP_OBO_URL)
    request.raise_for_status()

    for line in tqdm(request.text.splitlines(), unit=" lines"):
        line = line.rstrip("\n")
        value = " ".join(line.split(" ")[1:])
        if line.startswith("id: "):
            hpo_id = value
            hpo_id_to_record[hpo_id] = {
                'hpo_id': hpo_id,
                'is_category': False,
            }
        elif line.startswith("is_a: "):
            is_a = value.split(" ! ")[0]
            if is_a == "HP:0000118":
                hpo_id_to_record[hpo_id]['is_category'] = True
            hpo_id_to_record[hpo_id]['parent_id'] = is_a
        elif line.startswith("name: "):
            hpo_id_to_record[hpo_id]['name'] = value
        elif line.startswith("def: "):
            hpo_id_to_record[hpo_id]['definition'] = value
        elif line.startswith("comment: "):
            hpo_id_to_record[hpo_id]['comment'] = value

    return hpo_id_to_record


def parse_hpo_terms_arg(hpo_terms_arg, hpo_id_to_record):
    results = []
    skipped_counter = 0
    for hpo_terms in hpo_terms_arg:
        hpo_terms = hpo_terms.split(",")
        for hpo_term in hpo_terms:
            hpo_term = hpo_term.strip()
            if not hpo_term.startswith("HP:"):
                try:
                    hpo_term = f"HP:{int(hpo_term)}"
                except ValueError:
                    print(f"WARNING: Invalid HPO term: '{hpo_term}'. Skipping...")
                    skipped_counter += 1
                    continue

            if hpo_term not in hpo_id_to_record:
                print(f"WARNING: HPO term '{hpo_term}' not found in hp.obo. Skipping...")
                skipped_counter += 1
                continue

            results.append(hpo_term)

    if skipped_counter > 0:
        print(f"Skipped {skipped_counter} invalid HPO terms")

    return results


def get_category_id(hpo_id_to_record, hpo_id):
    """For a given hpo_id, get the hpo id of it's top-level category (eg. 'cardiovascular') and
    return it. If the hpo_id belongs to multiple top-level categories, return one of them.
    """

    if hpo_id == "HP:0000001":
        return None

    if 'parent_id' not in hpo_id_to_record[hpo_id]:
        return None

    while hpo_id_to_record[hpo_id]['parent_id'] != "HP:0000118":

        hpo_id = hpo_id_to_record[hpo_id]['parent_id']
        if hpo_id == "HP:0000001":
            return None
        if hpo_id not in hpo_id_to_record:
            raise ValueError("Unknown HPO id: %s" % hpo_id)

    return hpo_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", help="Print HPO term definitions")
    parser.add_argument("hpo_terms", nargs="+", help="comma- or space-separated list of HPO terms like HP:5200135")
    args = parser.parse_args()

    hpo_id_to_record = download_hpo_obo_file()

    # for each hpo id, find its top level category
    for hpo_id in hpo_id_to_record.keys():
        record = hpo_id_to_record[hpo_id]
        record["category_id"] = get_category_id(hpo_id_to_record, hpo_id)
        record["category"] = hpo_id_to_record[record["category_id"]]["name"] if record["category_id"] else ""

    print("Parsed %d HPO terms" % len(hpo_id_to_record))

    hpo_terms = parse_hpo_terms_arg(args.hpo_terms, hpo_id_to_record)
    print(f"{len(hpo_terms):,d} HPO terms:")

    category_field_width = max(len(hpo_id_to_record[hpo_term]["category"]) for hpo_term in hpo_terms)
    for hpo_term in sorted(hpo_terms, key=lambda hpo_term: (hpo_id_to_record[hpo_term]['category_id'], hpo_term)):
        record = hpo_id_to_record[hpo_term]

        if args.verbose:
            definition = record['definition'].split("[")[0].strip('" .').replace('\\"', '"')
            definition = f"    ({definition})"
        else:
            category = record["category"]
            category = f"{category:{category_field_width}s}"
            definition = ""

        print(f"{hpo_term}  {category} :    {record['name']}{definition}")


if __name__ == "__main__":
    main()