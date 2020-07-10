from functions import get_uniprot_entity, restAPI_get, get_fasta_path
from AA_properties import *

def get_feature(name_of_feature, data_dir, pdb_id, chain_id):
    if name_of_feature == 'PTM':
        return get_PTM(pdb_id, chain_id)
    elif name_of_feature == "hydropathy":
        return get_hydropathy(data_dir, pdb_id, chain_id)
    else:
        print("Error: unknown feature")
        #todo


def get_PTM(pdb_id, chain_id):
    try:
        entities = get_uniprot_entity(pdb_id, chain_id)
    except:
        print("error")
        #sys.exit(3)
        #todo

    feature_vals = []

    for entity in entities:
        uniprot_id = entity[0]
        entity_id = entity[1]
        unp_start = entity[2]
        unp_end = entity[3]
        start_res_num = entity[4]
        end_res_num = entity[5]

        # uniprot - jen MOD_RES type (podkategorie PTM)
        # url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types=MOD_RES"

        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=PTM"

        response = restAPI_get(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            type = feature["type"]
            if (type == "DISULFID"):  # disulfide bond - 'begin' and 'end' mean connected AAs
                rng = {feat_begin, feat_end}
            else:
                rng = range(feat_begin, feat_end + 1)
            for i in rng:
                res = i - unp_start  # mapping pdb residues to uniprot entry, counting from 0
                if res >= 0 and res < unp_end - unp_start + 1:
                    feature_vector[res] = 1

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

    return feature_vals

def get_hydropathy(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    with open(fasta_file, 'r') as file:
        seq = file.read()

    result = get_AA_scores(hydropathy_kyte_doolitle, seq)

    return result



