from Bio import SeqIO

from helper import get_uniprot_entity, restAPI_get_json, get_fasta_path
from AA_properties import *
import random

def get_feature(name_of_feature, data_dir, pdb_id, chain_id):
    if name_of_feature == 'unp_PTM':
        return get_PTM(data_dir, pdb_id, chain_id)
    elif name_of_feature == 'unp_glycosylation':
        return get_glycosylation(data_dir, pdb_id, chain_id)
    elif name_of_feature == 'unp_variants':
        return get_variants(data_dir, pdb_id, chain_id)
    elif name_of_feature == 'pdbkb_conservation':
        return get_pdbkb_conservation(data_dir, pdb_id, chain_id)
    #elif name_of_feature == 'unp_metal': #todo smazat, blbost
    #    return get_metal_binding(data_dir, pdb_id, chain_id)
    elif name_of_feature == "hydropathy":
        return get_hydropathy_kyte_doolitle(data_dir, pdb_id, chain_id)
    elif name_of_feature == "pKa_COOH":
        return get_pKa_COOH(data_dir, pdb_id, chain_id)
    elif name_of_feature == "pKa_NH3":
        return get_pKa_NH3(data_dir, pdb_id, chain_id)
    elif name_of_feature == "molecular_weight":
        return get_molecular_weight(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "random": #todo smazat
    #    return get_random(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "AA": #todo smazat
    #    return get_AAs(data_dir, pdb_id, chain_id)
    elif name_of_feature == "dynamine":
        return get_dynamine(data_dir, pdb_id, chain_id)
    else:
        print(f"Error: unknown feature {name_of_feature}.") #todo
        return

types_of_features = {
    "unp_PTM" : "discrete",
    "unp_glycosylation" : "discrete",
    "unp_variants" : "discrete",
    "hydropathy" : "continuous",
    "molecular_weight" : "continuous",
    "dynamine" : "continuous",
    "pKa_COOH" : "continuous",
    "pKa_NH3" : "continuous"
}

default_values = {
    "unp_PTM" : 0,
    "unp_glycosylation" : 0,
    "unp_variants" : 0,
    "hydropathy" : 0,
    "molecular_weight" : 110,
    "dynamine" : 5,
    "pKa_COOH" : 9.5,
    "pKa_NH3" : 2.2,
    "pdbkb_conservation" : 0
} #todo

def get_hydropathy_kyte_doolitle(data_dir, pdb_id, chain_id):
    return get_AA_properties(hydropathy_kyte_doolitle, data_dir, pdb_id, chain_id)

def get_molecular_weight(data_dir, pdb_id, chain_id):
    return get_AA_properties(molecular_weight, data_dir, pdb_id, chain_id)

def get_pKa_COOH(data_dir, pdb_id, chain_id):
    return get_AA_properties(pKa_COOH, data_dir, pdb_id, chain_id)

def get_pKa_NH3(data_dir, pdb_id, chain_id):
    return get_AA_properties(pKa_NH3, data_dir, pdb_id, chain_id)

def get_PTM(data_dir, pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)

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

        response = restAPI_get_json(url)

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

def get_glycosylation(data_dir, pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)

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

        response = restAPI_get_json(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            type = feature["type"]
            if (type != "CARBOHYD"):
                continue
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

def get_variants(data_dir, pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)

    feature_vals = []

    for entity in entities:
        uniprot_id = entity[0]
        entity_id = entity[1]
        unp_start = entity[2]
        unp_end = entity[3]
        start_res_num = entity[4]
        end_res_num = entity[5]

        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=VARIANTS"

        response = restAPI_get_json(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
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

def get_pdbkb_conservation(data_dir, pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)

    feature_vals = []

    for entity in entities:
        uniprot_id = entity[0]
        entity_id = entity[1]
        unp_start = entity[2]
        unp_end = entity[3]
        start_res_num = entity[4]
        end_res_num = entity[5]

        url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/sequence_conservation/{pdb_id}/{entity_id}"

        response = restAPI_get_json(url)

        feature_vector = [0] * (end_res_num - start_res_num + 1)  # including both start and end AAs

        for resi in response[pdb_id]["data"]:
            feat_begin = int(resi["start"])
            feat_end = int(resi["end"])
            if (feat_begin != feat_end):
                raise ValueError("Pdb KB conservation: start is not same as end!!") #todo only for debug

            feature_vector[feat_begin-1] = resi["conservation_score"]

    return feature_vector

def get_metal_binding(data_dir, pdb_id, chain_id):

    entities = get_uniprot_entity(pdb_id, chain_id)

    feature_vals = []

    for entity in entities:
        uniprot_id = entity[0]
        entity_id = entity[1]
        unp_start = entity[2]
        unp_end = entity[3]
        start_res_num = entity[4]
        end_res_num = entity[5]

        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types=METAL"

        response = restAPI_get_json(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
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

def get_AAs(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    with open(fasta_file, 'r') as file:
        seq = file.read()
    result = []
    for i in range(1, len(seq) + 1):
        AA = seq[i-1]
        result.append((i, AA))
    return result

def get_random(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    with open(fasta_file, 'r') as file:
        seq = file.read()
    result = []
    for i in range(1, len(seq) + 1):
        AA_score = random.uniform(1.0, 3.0)
        result.append((i, AA_score))
    return result

def get_AA_properties(scores_dict, data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    return get_AA_scores(scores_dict, seq)

def get_dynamine(data_dir, pdb_id, chain_id):
    import os
    from DynaMine.predictor import DynaMine

    result_dir = f"{data_dir}dynamine/"
    fasta_path = get_fasta_path(data_dir, pdb_id, chain_id)

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    dynamine = DynaMine(result_dir)
    if dynamine.predict(fasta_path, pdb_id, chain_id):
        print(f"DynaMine: {pdb_id} {chain_id} successfully processed.")
        result = []
        res_num = 1
        with open(f"{data_dir}/dynamine/{pdb_id}{chain_id}.txt") as f: #todo ukladat si rovnoujen ten soubor co potrebuju
            for line in f:
                if (line[0] != '*'):
                    val = float(line.split()[1])
                    result.append((res_num, val))
                    res_num += 1
        return result
    else:
        print(f"Error: DynaMine: {pdb_id} {chain_id}")


###### DEBUGGING ###############3
pdb_id = "1cbs"
chain_id = "A"
get_pdbkb_conservation(None, pdb_id, chain_id)