from Bio import SeqIO
from Bio.PDB import PDBParser, HSExposure, Selection, is_aa
from helper import *
from AA_properties import *
import random

def get_feature(name_of_feature, data_dir, pdb_id, chain_id):
    if name_of_feature == 'unp_PTM':
        return get_unp_PTM(pdb_id, chain_id)
    elif name_of_feature == 'unp_glycosylation':
        return get_unp_glycosylation(pdb_id, chain_id)
    elif name_of_feature == 'unp_lipidation':
        return get_unp_lipidation(pdb_id, chain_id)
    elif name_of_feature == 'unp_mod_res':
        return get_unp_mod_res(pdb_id, chain_id)
    elif name_of_feature == 'unp_variation':
        return get_unp_variation(pdb_id, chain_id)
    elif name_of_feature == 'unp_topology':
        return get_unp_topology(pdb_id, chain_id)
    elif name_of_feature == 'unp_compbias':
        return get_unp_compbias(pdb_id, chain_id)
    elif name_of_feature == 'pdbkb_sec_str':
        return get_pdbkb_sec_str(pdb_id, chain_id)
    #elif name_of_feature == 'unp_variants': #todo
    #    return get_variants(data_dir, pdb_id, chain_id)
    elif name_of_feature == 'pdbkb_conservation':
        return get_pdbkb_conservation(pdb_id, chain_id)
    elif name_of_feature == "aa":
        return get_aa(data_dir, pdb_id, chain_id)
    elif name_of_feature == "aa_hydropathy":
        return get_hydropathy_kyte_doolitle(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "aa_pKa_COOH":
    #    return get_pKa_COOH(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "aa_pKa_NH3":
    #    return get_pKa_NH3(data_dir, pdb_id, chain_id)
    elif name_of_feature == "aa_molecular_weight":
        return get_molecular_weight(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "random": #todo debug
    #    return get_random(data_dir, pdb_id, chain_id)
    #elif name_of_feature == "AA": #todo debug
    #    return get_AAs(data_dir, pdb_id, chain_id)
    elif name_of_feature == "dynamine":
        return get_dynamine(data_dir, pdb_id, chain_id)
    elif name_of_feature == "mobiDB":
        return get_mobiDB_lite(pdb_id, chain_id)
    elif name_of_feature == "HSE":
        return get_HSE(data_dir, pdb_id, chain_id)
    elif name_of_feature == "HSE_down":
        return get_HSE_down(data_dir, pdb_id, chain_id)
    elif name_of_feature == "exposureCN":
        return get_exposureCN(data_dir, pdb_id, chain_id)
    elif name_of_feature == "b_factor":
        return get_b_factor(data_dir, pdb_id, chain_id)
    else:
        raise ValueError(f"Unknown feature {name_of_feature}.") #todo
        return

 #todo aby tohle nebylo potreba
types_of_features = {
    "unp_PTM" : "binary",
    "unp_glycosylation" : "binary",
    "unp_lipidation" : "binary",
    "unp_mod_res" : "binary",
    "unp_variation" : "binary",
    "unp_topology" : "binary",
    "unp_compbias" : "binary",
    "pdbkb_sec_str" : "categorical",
    "aa_hydropathy" : "continuous",
    "aa_molecular_weight" : "continuous",
    "dynamine" : "continuous",
    "pdbkb_conservation" : "categorical",
    "mobiDB" : "continuous",
    "aa" : "categorical",
    "HSE" : "continuous",
    "HSE_down" : "continuous",
    "exposureCN" : "continuous",
    "b_factor" : "continuous"
}

default_values = {
    "unp_PTM" : 0,
    "unp_glycosylation" : 0,
    "unp_lipidation" : 0,
    "unp_mod_res" : 0,
    "unp_variants" : 0,
    "hydropathy" : 0,
    "molecular_weight" : 110,
    "dynamine" : 0.8,
    "pKa_COOH" : 9.5,
    "pKa_NH3" : 2.2,
    "pdbkb_conservation" : 0
} #todo




def get_uniprot_type(type, pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)
    feature_vals = []
    for entity in entities:
        uniprot_id = entity[0]
        unp_start = entity[1]
        unp_end = entity[2]
        start_res_num = entity[3]
        end_res_num = entity[4]
        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types={type}"
        response = restAPI_get_json(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs
        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            for i in range(feat_begin, feat_end + 1):
                res = i - unp_start  # mapping pdb residues to uniprot entry, counting from 0
                if res >= 0 and res < unp_end - unp_start + 1:
                    feature_vector[res] = 1

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

    return feature_vals

def get_unp_glycosylation(pdb_id, chain_id):
    return get_uniprot_type("CARBOHYD", pdb_id, chain_id)

def get_unp_lipidation(pdb_id, chain_id):
    return get_uniprot_type("LIPID", pdb_id, chain_id)

def get_unp_mod_res(pdb_id, chain_id):
    return get_uniprot_type("MOD_RES", pdb_id, chain_id)

def get_unp_PTM(pdb_id, chain_id):
    return get_uniprot_type("CARBOHYD%2CMOD_RES%2CLIPID", pdb_id, chain_id) #glycosylation+lipidation+mod_res

def get_unp_topology(pdb_id, chain_id):
    return get_uniprot_type("TRANSMEM%2CINTRAMEM", pdb_id, chain_id) #transmembrane + intramembrane

def get_unp_compbias(pdb_id, chain_id):
    return get_uniprot_type("COMPBIAS", pdb_id, chain_id)

def get_unp_variation(pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)
    feature_vals = []
    for entity in entities:
        uniprot_id = entity[0]
        unp_start = entity[1]
        unp_end = entity[2]
        start_res_num = entity[3]
        end_res_num = entity[4]

        url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
        response = restAPI_get_json(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            if (feat_begin != feat_end):
                raise ValueError(f"ERROR: {uniprot_id}: url = {url} - feat_begin != feat_end") #todo only for debugging
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

def get_AA_properties(scores_dict, data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    return get_AA_scores(scores_dict, seq)

def get_hydropathy_kyte_doolitle(data_dir, pdb_id, chain_id):
    return get_AA_properties(hydropathy_kyte_doolitle, data_dir, pdb_id, chain_id)

def get_molecular_weight(data_dir, pdb_id, chain_id):
    return get_AA_properties(molecular_weight, data_dir, pdb_id, chain_id)

def get_pKa_COOH(data_dir, pdb_id, chain_id):
    return get_AA_properties(pKa_COOH, data_dir, pdb_id, chain_id)

def get_pKa_NH3(data_dir, pdb_id, chain_id):
    return get_AA_properties(pKa_NH3, data_dir, pdb_id, chain_id)



def get_dynamine(data_dir, pdb_id, chain_id): #todo cely upravit, zkontrolovat
    import os
    from DynaMine.predictor import DynaMine

    result_dir = os.path.join(data_dir, "dynamine")
    fasta_path = get_fasta_path(data_dir, pdb_id, chain_id)

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    dynamine = DynaMine(result_dir)
    if dynamine.predict(fasta_path, pdb_id, chain_id):
        #print(f"DynaMine: {pdb_id} {chain_id} successfully processed.")
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
        raise ValueError(f"Error: DynaMine: {pdb_id} {chain_id}")

def get_pdbkb_conservation(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/sequence_conservation/{pdb_id}/{entity_id}"
    response = restAPI_get_json(url)

    feature_vals = []
    for resi in response[pdb_id]["data"]:
        feat_begin = int(resi["start"])
        feat_end = int(resi["end"])
        if (feat_begin != feat_end):
            raise ValueError("Pdb KB conservation: start is not same as end!!") #todo only for debug
        feature_vals.append((feat_begin, resi["conservation_score"]))

    return feature_vals

def get_pdbkb_sec_str(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/secondary_structure/{pdb_id}/{entity_id}"
    response = restAPI_get_json(url)

    feature_vector = ['N'] * int(response[pdb_id]["length"])
    for rec in response[pdb_id]["data"]:
        dataType = rec["dataType"] #Helix, Strand, MobiDB
        if (dataType == "Helix"):
            val = 'H'
        elif (dataType == "Strand"):
            val = 'S'
        elif (dataType == "MobiDB"): #todo
            val = 'M'
        else:
            raise ValueError(f"PDBe KB - secondary structure: {pdb_id} {chain_id} : unknown dataType {dataType}")
        for resi in rec["residues"]:
            feat_begin = int(resi["startIndex"])
            feat_end = int(resi["endIndex"])
            for i in range(feat_begin, feat_end + 1):
                feature_vector[i-1] = val

    feature_vals = []
    for res_num in range(1, int(response[pdb_id]["length"] + 1)):
        feature_vals.append((res_num, feature_vector[res_num-1]))

    return feature_vals


def get_aa(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    feature_vals = []
    for i in range(1, len(seq) + 1):
        AA = seq[i-1]
        feature_vals.append((i, AA))
    return feature_vals

'''def get_random(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    with open(fasta_file, 'r') as file:
        seq = file.read()
    result = []
    for i in range(1, len(seq) + 1):
        AA_score = random.uniform(1.0, 3.0)
        result.append((i, AA_score))
    return result
'''

def get_mobiDB_lite(pdb_id, chain_id):
    entities = get_uniprot_entity(pdb_id, chain_id)
    feature_vals = []
    for entity in entities:
        uniprot_id = entity[0]
        unp_start = entity[1]
        unp_end = entity[2]
        start_res_num = entity[3]
        end_res_num = entity[4]

        url = f"https://mobidb.bio.unipd.it/ws/{uniprot_id}/consensus"
        response = restAPI_get_json(url)

        #feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs #todo default value misto 0

        pred = response["mobidb_consensus"]["disorder"]["predictors"]
        for p in pred:
            if (p["method"] == "mobidb_lite"):
                vals = p["scores"]
                break

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, vals[unp_start + i - 1]))
            i += 1

    return feature_vals

def get_HSE(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]
    HSExposure.HSExposureCB(model)
    feature_vals = []
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    for r in model.get_residues():
        if not (r.id[0].isspace()): #is HETATM
            continue
        try:
            hse = abs(int(r.xtra["EXP_HSE_B_U"]) - int(r.xtra["EXP_HSE_B_D"]))
            auth_res_num = getFullAuthorResNum(r.id)
            #auth_res_num = str(r.id[1])
            #if (r.id[2] != ' '):
            #   auth_res_num += str(r.id[2])  # insertion code
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, hse))
        except:
            print(f"Error HSE: {pdb_id} {chain_id}: residue {r}") #todo

    return feature_vals

def get_HSE_down(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]
    HSExposure.HSExposureCB(model)
    feature_vals = []
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    for r in model.get_residues():
        if not (r.id[0].isspace()):  # is HETATM
            continue
        try:
            hse = int(r.xtra["EXP_HSE_B_D"])
            auth_res_num = getFullAuthorResNum(r.id)
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, hse))
        except:
            pass
            #print(f"Error HSE: {pdb_id} {chain_id}: residue {r}") #todo

    return feature_vals

def get_exposureCN(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]
    HSExposure.ExposureCN(model)
    feature_vals = []
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    for r in model.get_residues():
        if not (r.id[0].isspace()):  # is HETATM
            continue
        try:
            cn = int(r.xtra["EXP_CN"])
            auth_res_num = getFullAuthorResNum(r.id)
            #auth_res_num = str(r.id[1])
            #if (r.id[2] != ' '):
            #    auth_res_num += str(r.id[2])  # insertion code
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, cn))
        except:
            print(f"Error CN: {pdb_id} {chain_id}: residue {r}") #todo

    return feature_vals

def get_b_factor(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, pdb_file)  #todo udelat funkci get_chain?
    chain = structure[0][chain_id]
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    feature_vals = []
    for residue in chain.get_residues():
        if not isPartOfChain(residue, mappings):
            continue
        sum = 0
        for atom in residue.child_list:
            sum += atom.bfactor
        b_factor = sum / len(residue.child_list)
        auth_res_num = getFullAuthorResNum(residue.id)
        pdbe_res_num = mappings[auth_res_num]
        feature_vals.append((pdbe_res_num, b_factor))
    return feature_vals


###DEBGUG###
#data_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_test"
#pdb_id = "1qhi"
#chain_id = "A"
#print(get_b_factor(data_dir, pdb_id, chain_id))
