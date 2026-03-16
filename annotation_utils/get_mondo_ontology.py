from annotation_utils.cache_utils import cache_json
import lxml.etree
import os
import requests

MONDO_RARE_OWL_URL = "https://purl.obolibrary.org/obo/mondo/subsets/mondo-rare.owl"
MONDO_OBO_URL = "https://purl.obolibrary.org/obo/mondo.obo"

@cache_json
def get_mondo_rare_disease_terms(mondo_rare_owl_file_path):
    """Get rare disease terms from the mondo-rare.owl file which is a separate download from the main mondo.obo file

    Note: This list includes cancer and rare disease terms that are not marked as 'rare' in the mondo.obo file.

    Args:
        mondo_rare_owl_file_path: Path to the mondo-rare.owl file
    Returns:
        List of mondo rare disease terms (eg. "MONDO:0958331")
    """
    mondo_rare_disease_terms = []
    
    for _, class_element in lxml.etree.iterparse(mondo_rare_owl_file_path, events=('end',), tag='{*}Class'):
        for key, value in class_element.items():
            if key.endswith("about"):
                mondo_term_url = value
                mondo_term = os.path.basename(mondo_term_url)
                mondo_rare_disease_terms.append(mondo_term.replace("_", ":")) 
                break

        # Clear element to free memory
        class_element.clear()
    
    return mondo_rare_disease_terms  # Example: MONDO:0958331

@cache_json
def download_mondo_obo_file():
    """    
    Parse the mondo.obo file and return a dictionary that maps mondo ids to a record containing the following fields:
    - mondo_id
    - name
    - definition
    - xrefs
    - is_rare
    - parent_id

    Example:
    {
        "MONDO:8000011": {
            "mondo_id": "MONDO:8000011",
            "name": "visceral neuropathy, familial, 1, autosomal recessive",
            "definition": "A form of chronic intestinal pseudoobstruction caused by a developmental failure of the enteric neurons to differentiate or migrate properly and manifests as a bowel obstruction.",
            "xrefs": ["DOID:0080679", "GARD:3928", "MEDGEN:340946", "MESH:C537394", "OMIM:243180", "Orphanet:99811", "UMLS:C1855733"],
            "is_rare": True,
            "parent_id": "MONDO:0000858"
        }
    }

    Returns:
        Dictionary that maps mondo ids to a record containing the above fields

    Example of a record in the mondo.obo file:

    [Term]
    id: MONDO:8000011
    name: visceral neuropathy, familial, 1, autosomal recessive
    def: "A form of chronic intestinal pseudoobstruction caused by a developmental failure of the enteric neurons to differentiate or migrate properly and manifests as a bowel obstruction." [Orphanet:99811]
    subset: gard_rare {source="GARD:3928", source="MONDO:GARD"}
    subset: nord_rare {source="MONDO:NORD"}
    subset: ordo_subtype_of_a_disorder {source="Orphanet:99811"}
    subset: rare
    synonym: "Argyrophil myenteric plexus deficiency of" RELATED [GARD:0003969]
    synonym: "Argyrophil myenteric plexus, deficiency of" RELATED []
    synonym: "intestinal pseudoobstruction due to neuronal disease" RELATED []
    synonym: "neuronal intestinal dysplasia, type a" RELATED []
    synonym: "NID A" RELATED []
    synonym: "pseudoobstruction chronic idiopathic intestinal neuronal type" RELATED [GARD:0003969]
    synonym: "pseudoobstruction, chronic idiopathic intestinal, neuronal type" RELATED []
    synonym: "visceral neuropathy familial" RELATED [GARD:0003969]
    synonym: "visceral neuropathy, familial, autosomal recessive" RELATED []
    xref: DOID:0080679 {source="MONDO:equivalentTo"}
    xref: GARD:3928 {source="MONDO:GARD"}
    xref: MEDGEN:340946 {source="MONDO:equivalentTo", source="MONDO:MEDGEN"}
    xref: MESH:C537394 {source="Orphanet:99811", source="MONDO:equivalentTo", source="Orphanet:99811/e"}
    xref: OMIM:243180 {source="MONDO:equivalentTo"}
    xref: Orphanet:99811 {source="MONDO:equivalentTo"}
    xref: UMLS:C1855733 {source="MEDGEN:340946", source="MONDO:equivalentTo", source="MONDO:MEDGEN"}
    is_a: MONDO:0000858 {source="DC-OMIM:243180", source="DOID:0080679"} ! neuronal intestinal dysplasia
    is_a: MONDO:0017574 {source="Orphanet:99811"} ! chronic intestinal pseudoobstruction
    is_a: MONDO:0023961 {source="OMIM:243180"} ! visceral neuropathy, familial
    relationship: curated_content_resource https://www.malacards.org/card/neuronal_intestinal_dysplasia_type_a {source="MONDO:MalaCards"}
    relationship: curated_content_resource https://www.malacards.org/card/visceral_neuropathy_familial_1_autosomal_recessive {source="MONDO:MalaCards"}
    relationship: has_material_basis_in_germline_mutation_in http://identifiers.org/hgnc/3431 {source="OMIM:243180"} ! ERBB3
    property_value: skos:exactMatch DOID:0080679
    property_value: skos:exactMatch http://identifiers.org/medgen/340946
    property_value: skos:exactMatch http://identifiers.org/mesh/C537394
    property_value: skos:exactMatch http://linkedlifedata.com/resource/umls/id/C1855733
    property_value: skos:exactMatch https://omim.org/entry/243180
    property_value: skos:exactMatch Orphanet:99811
    property_value: terms:creator https://orcid.org/0000-0001-5208-3432
    """

    # top level parent of disease ontology: MONDO:0700096   (Human Disease). It has 42 child terms (eg. MONDO:0002051 integumentary system disorder)

    print("Downloading mondo.obo data")

    mondo_id = None
    mondo_id_to_record = {}
    request = requests.get(MONDO_OBO_URL)
    request.raise_for_status()

    is_inside_term_record = False
    for line in request.text.splitlines():
        line = line.rstrip("\n")
        if line.strip() == "":
            is_inside_term_record = False
            continue
        if line.startswith("[Term]"):
            mondo_id = None
            is_inside_term_record = True
            continue
        if not is_inside_term_record:
            continue

        value = " ".join(line.split(" ")[1:])
        value = value.split(" {")[0]

        if line.startswith("id: "):
            if not value.startswith("MONDO:"):
                is_inside_term_record = False
                continue

            mondo_id = value
            mondo_id_to_record[mondo_id] = {
                'mondo_id': mondo_id,
                'parent_id': None,
                'is_category': False,
                'xrefs': [],
            }
        elif line.startswith("is_a: "):
            is_a = value.split(" ! ")[0]
            if not is_a.startswith("MONDO:"):
                continue
            if is_a == "MONDO:0700096":
                mondo_id_to_record[mondo_id]['is_category'] = True
            mondo_id_to_record[mondo_id]['parent_id'] = is_a
        elif line.startswith("name: "):
            mondo_id_to_record[mondo_id]['name'] = value
        elif line.startswith("subset: ") and value == "rare":
            mondo_id_to_record[mondo_id]['is_rare'] = True
        elif line.startswith("def: "):
            mondo_id_to_record[mondo_id]['definition'] = value
        elif line.startswith("xref: "):
            mondo_id_to_record[mondo_id]['xrefs'].append(value)

    return mondo_id_to_record



def get_category_id(mondo_id_to_record, mondo_id):
    """For a given mondo_id, get the mondo id of it's top-level category (eg. 'cardiovascular') and
    return it. If the mondo_id belongs to multiple top-level categories, return one of them.
    """

    if mondo_id == "MONDO:0700096":
        return None

    if 'parent_id' not in mondo_id_to_record[mondo_id]:
        return None

    while mondo_id_to_record[mondo_id].get('parent_id') != "MONDO:0700096":

        if not mondo_id_to_record[mondo_id].get('parent_id'):
            #print(f"WARNING: No parent_id found for {mondo_id_to_record[mondo_id]}")
            return None

        mondo_id = mondo_id_to_record[mondo_id].get('parent_id')
        if mondo_id == "MONDO:0700096":
            return None
        if mondo_id not in mondo_id_to_record:
            raise ValueError("Unknown MONDO id: %s" % mondo_id)

    return mondo_id



def get_mondo_ontology():
    """
    # Ignore the mondo-rare.owl file because it is redundant with the mondo.obo file's "subset: rare" annotation
    mondo_rare_owl_file_path = "mondo-rare.owl"
    if not os.path.exists(mondo_rare_owl_file_path):
        print("Downloading mondo-rare.owl")
        os.system(f"wget {MONDO_RARE_OWL_URL}")    

    print("Parsing mondo-rare.owl")
    mondo_rare_disease_terms = get_mondo_rare_disease_terms(mondo_rare_owl_file_path)
    print(f"Parsed {len(mondo_rare_disease_terms):,d} mondo rare disease terms from {mondo_rare_owl_file_path}")
    assert len(mondo_rare_disease_terms) == len(set(mondo_rare_disease_terms)), "Duplicate mondo terms found in mondo-rare.owl"
    """

    # parse the OBO file with all terms and subset it to the rare disease terms
    mondo_term_lookup = download_mondo_obo_file()

    mondo_rare_disease_term_lookup = {
        mondo_id: record for mondo_id, record in mondo_term_lookup.items() if record.get('is_rare')
    }

    # for each mondo term, record its top level category
    for mondo_id in mondo_rare_disease_term_lookup.keys():
        record = mondo_term_lookup[mondo_id]
        record["category_id"] = get_category_id(mondo_term_lookup, mondo_id)
        record["category"] = mondo_term_lookup[record["category_id"]]["name"] if record["category_id"] else ""

    print(f"Parsed {len(mondo_term_lookup):,d} mondo terms from mondo.obo, of which {len(mondo_rare_disease_term_lookup):,d} are rare disease terms")


    #print(f"Parsed {len(mondo_lookup):,d} mondo terms from mondo.obo, {len(mondo_rare_disease_term_lookup):,d} of which are rare disease terms, and {len(set(mondo_rare_disease_terms) & set(mondo_rare_disease_term_lookup.keys())):,d} of which are in the mondo-rare.owl file")
    #print(set(mondo_rare_disease_terms) - set(mondo_rare_disease_term_lookup.keys()))

    return mondo_rare_disease_term_lookup

if __name__ == "__main__":
    get_mondo_ontology()
