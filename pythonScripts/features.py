from Bio import SeqIO
from Bio.PDB import PDBParser, HSExposure, Selection, is_aa
from helper import *
from AA_properties import *
import os
import random

def get_feature(name_of_feature, data_dir, pdb_id, chain_id):
    if name_of_feature == 'unp_PTM':
        return get_uniprot_binary_type("CARBOHYD%2CMOD_RES%2CLIPID", pdb_id, chain_id) #glycosylation+lipidation+mod_res
    elif name_of_feature == 'unp_glycosylation':
        return get_uniprot_binary_type("CARBOHYD", pdb_id, chain_id)
    elif name_of_feature == 'unp_lipidation':
        return get_uniprot_binary_type("LIPID", pdb_id, chain_id)
    elif name_of_feature == 'unp_mod_res':
        return get_uniprot_binary_type("MOD_RES", pdb_id, chain_id)
    elif name_of_feature == 'unp_variation':
        return get_unp_variation(pdb_id, chain_id)
    elif name_of_feature == 'unp_topology':
        return get_uniprot_binary_type("TRANSMEM%2CINTRAMEM", pdb_id, chain_id) #transmembrane + intramembrane
    elif name_of_feature == 'unp_sec_str':
        return get_uniprot_sec_str(pdb_id, chain_id)
    elif name_of_feature == 'unp_non_standard':
        return get_uniprot_binary_type("NON_STD", pdb_id, chain_id)
    elif name_of_feature == 'unp_natural_variant':
        return get_uniprot_binary_type("VARIANT", pdb_id, chain_id)
    elif name_of_feature == 'unp_compbias':
        return get_uniprot_binary_type("COMPBIAS", pdb_id, chain_id)
    elif name_of_feature == 'pdbekb_sec_str':
        return get_pdbkb_sec_str(pdb_id, chain_id)
    elif name_of_feature == 'pdbekb_conservation':
        return get_pdbkb_conservation(pdb_id, chain_id)
    elif name_of_feature == "conservation":
        return get_conservation(data_dir, pdb_id, chain_id)
    elif name_of_feature == "aa":
        return get_aa(data_dir, pdb_id, chain_id)
    elif name_of_feature == "aa_pairs":
        return get_aa_pairs(data_dir, pdb_id, chain_id)
    elif name_of_feature == "hydropathy":
        return get_AA_properties(hydropathy_kyte_doolitle, data_dir, pdb_id, chain_id)
    elif name_of_feature == "polarity":
        return get_AA_properties(polarity, data_dir, pdb_id, chain_id)
    elif name_of_feature == "polarity_binary":
        return get_AA_properties(polarity_binary, data_dir, pdb_id, chain_id)
    elif name_of_feature == "charged":
        return get_AA_properties(charged, data_dir, pdb_id, chain_id)
    elif name_of_feature == "aromaticity":
        return get_AA_properties(aromaticity, data_dir, pdb_id, chain_id)
    elif name_of_feature == "mol_weight":
        return get_AA_properties(molecular_weight, data_dir, pdb_id, chain_id)
    elif name_of_feature == "H_bond_atoms":
        return get_AA_properties(H_bond_atoms, data_dir, pdb_id, chain_id)
    elif name_of_feature == "dynamine_website":
        return get_dynamine(data_dir, pdb_id, chain_id)
    elif name_of_feature == "dynamine_funPDBe":
        return get_funPDBe("dynamine", "backbone", pdb_id, chain_id)
    elif name_of_feature == "efoldmine_funPDBe":
        return get_funPDBe("dynamine", "efoldmine", pdb_id, chain_id)
    elif name_of_feature == "mobiDB":
        return get_mobiDB_lite(pdb_id, chain_id)
    elif name_of_feature == "HSE":
        return get_HSE(data_dir, pdb_id, chain_id)
    elif name_of_feature == "HSE_down":
        return get_HSE_down(data_dir, pdb_id, chain_id)
    elif name_of_feature == "exposureCN":
        return get_exposureCN(data_dir, pdb_id, chain_id)
    elif name_of_feature == "bfactor":
        return get_bfactor(data_dir, pdb_id, chain_id)
    elif name_of_feature == "bfactor_Calpha":
        return get_bfactor_Calpha(data_dir, pdb_id, chain_id)
    elif name_of_feature == "depth":
        return get_funPDBe("depth", "complex_residue_depth", pdb_id, chain_id)
    elif name_of_feature == "phi_angle":
        return get_phi_angle(pdb_id, chain_id)
    elif name_of_feature == "psi_angle":
        return get_psi_angle(pdb_id, chain_id)
    elif name_of_feature == "cis_peptide":
        return get_cis_peptide(pdb_id, chain_id)
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
    "unp_non_standard" : "binary",
    "unp_sec_str" : "categorical",
    "unp_natural_variant" : "binary",
    "pdbekb_sec_str" : "categorical",
    "hydropathy" : "continuous",
    "mol_weight" : "continuous",
    "dynamine_website" : "continuous",
    "dynamine_funPDBe" : "continuous",
    "efoldmine_funPDBe" : "continuous",
    "pdbekb_conservation" : "categorical", #todo hodnoty 0-4 ? -> zkusit continuous
    "mobiDB" : "continuous",
    "aa" : "categorical",
    "aa_pairs" : "categorical",
    "HSE" : "continuous",
    "HSE_down" : "continuous",
    "exposureCN" : "continuous",
    "bfactor" : "continuous",
    "bfactor_Calpha" : "continuous",
    "polarity" : "categorical",
    "polarity_binary" : "binary",
    "charged" : "binary",
    "aromaticity" : "binary",
    "H_bond_atoms" : "categorical",
    "depth" : "continuous",
    "phi_angle" : "continuous",
    "psi_angle" : "continuous",
    "cis_peptide" : "binary",
    "conservation" : "continuous"
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




def get_uniprot_binary_type(type, pdb_id, chain_id):
    segments = get_uniprot_segments(pdb_id, chain_id)
    feature_vals = []
    for seg in segments:
        uniprot_id = seg[0]
        segment_begin = seg[1]
        segment_end = seg[2]
        start_res_num = seg[3]
        end_res_num = seg[4]
        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types={type}"
        response = restAPI_get_json(url)

        feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs
        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            for i in range(max(feat_begin, segment_begin), min(feat_end, segment_end) + 1): #take only parts of features that overlap with uniprot segment
                res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                feature_vector[res] = 1

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

    return feature_vals

def get_uniprot_sec_str(pdb_id, chain_id):
    segments = get_uniprot_segments(pdb_id, chain_id)
    feature_vals = []
    for seg in segments:
        uniprot_id = seg[0]
        segment_begin = seg[1]
        segment_end = seg[2]
        start_res_num = seg[3]
        end_res_num = seg[4]
        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=STRUCTURAL"
        response = restAPI_get_json(url)

        feature_vector = ['X'] * (segment_end - segment_begin + 1)  # including both start and end AAs
        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            type = feature["type"]
            if type == "HELIX":
                val = 'H'
            elif type == "STRAND":
                val = 'S'
            elif type == "TURN":
                val = 'T'
            else:
                raise ValueError(f"Unknown feature type '{type}'")
            for i in range(max(feat_begin, segment_begin), min(feat_end,segment_end) + 1):
                res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                feature_vector[res] = val

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

    return feature_vals

def get_unp_variation(pdb_id, chain_id):
    segments = get_uniprot_segments(pdb_id, chain_id)
    feature_vals = []
    for seg in segments:
        uniprot_id = seg[0]
        segment_begin = seg[1]
        segment_end = seg[2]
        start_res_num = seg[3]
        end_res_num = seg[4]

        url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
        response = restAPI_get_json(url)

        feature_vector = [0] * (segment_end - segment_begin + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            if (feat_begin != feat_end):
                raise ValueError(f"ERROR: {uniprot_id}: url = {url} - feat_begin != feat_end") #todo only for debugging
            rng = range(feat_begin, feat_end + 1)
            for i in rng:
                res = i - segment_begin  # mapping pdb residues to uniprot entry, counting from 0
                if res >= 0 and res < segment_end - segment_begin + 1:
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

def get_phi_angle(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
    response = restAPI_get_json(url)
    molecules = response[pdb_id]["molecules"]
    count = 0
    feature_vals = []
    for molecule in molecules:
        if molecule["entity_id"] == entity_id:
            for chain in molecule["chains"]:
                if (chain["chain_id"] == chain_id):
                    for resi in chain["models"][0]["residues"]: #todo: muzu brat automaticky prvni model?
                        val = resi["phi"]
                        if (val != None):
                            feature_vals.append((resi["residue_number"], val))
                    count += 1
    if count != 1: #todo smazat, jen debug
        raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")
        return

    return feature_vals

def get_psi_angle(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
    response = restAPI_get_json(url)
    molecules = response[pdb_id]["molecules"]
    count = 0
    feature_vals = []
    for molecule in molecules:
        if molecule["entity_id"] == entity_id:
            for chain in molecule["chains"]:
                if (chain["chain_id"] == chain_id):
                    for resi in chain["models"][0]["residues"]: #todo: muzu brat automaticky prvni model?
                        val = resi["psi"]
                        if (val != None):
                            feature_vals.append((resi["residue_number"], val))
                    count += 1
    if count != 1: #todo smazat, jen debug
        raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")
        return

    return feature_vals

def get_cis_peptide(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/api/validation/rama_sidechain_listing/entry/{pdb_id}"
    response = restAPI_get_json(url)
    molecules = response[pdb_id]["molecules"]
    feature_vals = []
    for molecule in molecules:
        if molecule["entity_id"] == entity_id:
            for chain in molecule["chains"]:
                if (chain["chain_id"] == chain_id):
                    for resi in chain["models"][0]["residues"]: #todo: muzu brat automaticky prvni model?
                        val = resi["cis_peptide"]
                        if (val == 'Y'):
                            feature_vals.append((resi["residue_number"], 1))
                        else:
                            feature_vals.append((resi["residue_number"], 0))
    return feature_vals

def get_pdbkb_sec_str(pdb_id, chain_id):
    entity_id = get_entity_id(pdb_id, chain_id)
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/secondary_structure/{pdb_id}/{entity_id}"
    response = restAPI_get_json(url)

    feature_vector = ['X'] * int(response[pdb_id]["length"])
    for rec in response[pdb_id]["data"]:
        dataType = rec["dataType"] #Helix, Strand, MobiDB
        if (dataType == "Helix"):
            val = 'H'
        elif (dataType == "Strand"):
            val = 'S'
        elif (dataType == "MobiDB"): #skip this type
            continue
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

def get_aa_pairs(data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    feature_vals = []
    for i in range(1, len(seq)):
        aa_1 = seq[i-1]
        aa_2 = seq[i]
        feature_vals.append((i, aa_1 + aa_2))
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
    entities = get_uniprot_segments(pdb_id, chain_id)
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
    model = structure[0] #todo muzu brat automaticky prvni model?
    HSExposure.HSExposureCB(model)
    feature_vals = []
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    for r in model.get_residues():
        if not (r.id[0].isspace()): #is HETATM
            continue
        try:
            hse = abs(int(r.xtra["EXP_HSE_B_U"]) - int(r.xtra["EXP_HSE_B_D"]))
            auth_res_num = getFullAuthorResNum(r.id)
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, hse))
        except:
            pass #some residues are incomplete in PDB file (i.e. CB atom is missing) and no HSE value is returned

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
            pass #some residues are incomplete in PDB file (i.e. CB atom is missing) and no HSE value is returned

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
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, cn))
        except:
            raise ValueError(f"Error CN: {pdb_id} {chain_id}: residue {r}") #todo jestli se to pro cely dataset nestane, tak asi muzu smazat try-except

    return feature_vals

def get_bfactor(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)
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
        b_factor = round(sum / len(residue.child_list), 3)
        auth_res_num = getFullAuthorResNum(residue.id)
        pdbe_res_num = mappings[auth_res_num]
        feature_vals.append((pdbe_res_num, b_factor))
    return feature_vals

def get_bfactor_Calpha(data_dir, pdb_id, chain_id):
    pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
    parser = PDBParser(PERMISSIVE=0, QUIET=1)
    structure = parser.get_structure(pdb_id, pdb_file)  #todo udelat funkci get_chain?
    chain = structure[0][chain_id]
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
    feature_vals = []
    for residue in chain.get_residues():
        if not isPartOfChain(residue, mappings):
            continue
        bfactor = None
        for atom in residue.child_list:
            if atom.id == 'CA':
                if (bfactor != None):
                    raise ValueError(f"Error: more C alpha atoms in {pdb_id} {chain_id} residue {residue.id}") #todo only debug
                bfactor = atom.bfactor
                #todo zrychlit, kdyz najdu C alpha -> break
        if bfactor == None:
            print(f"Error: no C alpha in {pdb_id} {chain_id} residue {residue.id}") # todo only debug
            continue
        auth_res_num = getFullAuthorResNum(residue.id)
        pdbe_res_num = mappings[auth_res_num]
        feature_vals.append((pdbe_res_num, bfactor))
    return feature_vals

def get_funPDBe(resource, label, pdb_id, chain_id):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/funpdbe_annotation/{resource}/{pdb_id}"
    feature_vals = {}
    response = restAPI_get_json(url)
    bug = False #todo odstranit, az opravi bug v pdbekb
    for rec in response[pdb_id][0]["annotations"]:
        if rec["label"] == label:
            residues = rec["site_residues"]
            for res in residues:
                if (res["chain_id"] != chain_id):
                    continue
                score = res["raw_score"]
                res_num = res["residue_number"]

                if res_num in feature_vals:
                    bug = True
                feature_vals[res_num] = score

    if bug:
        raise ValueError(f"ERROR: {pdb_id} {chain_id}: {resource}: more values for some residues.")

    feature_vals_list = [(k, v) for k, v in feature_vals.items()]
    return feature_vals_list

def get_conservation(data_dir, pdb_id, chain_id):
    filepath = os.path.join(data_dir, "conservation", pdb_id + chain_id + ".json")
    feature_vals = []
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id) #todo smazat
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0] #todo smazat, debug
    with open(filepath) as json_file:
        data = json.load(json_file)
        scores = data["conservation"]
        if (len(seq) != len(scores)):
            raise ValueError(f"ERROR!!!!! {pdb_id} {chain_id} length of conservation scores not same as fasta seq") #todo smazat
        for i in range(0, len(scores)):
            score = scores[i]
            if (score < 0):
                score = 0
            feature_vals.append((i+1, score))
    return feature_vals


##DEBUG###
#data_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_chen11"
#pdb_id = "2hyr"
#chain_id = "A"
#print(get_feature("mobiDB", data_dir, pdb_id, chain_id))
#print()