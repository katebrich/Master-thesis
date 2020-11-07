from __future__ import print_function

import os
import urllib.request
import json
import sys
import xmltodict
import numpy as np
from urllib.error import HTTPError

#returns uniprotID, entityID, start, end
from Bio.PDB import is_aa


def get_uniprot_entity(pdb_id, chain_id): #todo cache
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'
    parsedResponse = restAPI_get_json(url)
    uniprotRecords = parsedResponse[f"{pdb_id}"]["UniProt"]

    #find entity with given chain_id
    entities = [] #todo dictionary?
    for record in uniprotRecords:
        uniprot_ID = record
        mappings = uniprotRecords[record]["mappings"]
        for entity in mappings:
            if entity["chain_id"] == chain_id:
                start_res_num = entity["start"]["residue_number"]
                end_res_num = entity["end"]["residue_number"]
                unp_start =  entity["unp_start"]
                unp_end =  entity["unp_end"]
                if (unp_end - unp_start != end_res_num - start_res_num):
                    raise ValueError(
                        f"Incorrect mapping: {pdb_id} {chain_id}, {uniprot_ID}: unp {unp_start}-{unp_end}, pdbe {start_res_num}-{end_res_num}")  # todo debug
                entities.append((uniprot_ID,unp_start, unp_end, start_res_num, end_res_num))

    if len(entities) < 1:
        sys.stderr.write(f"Error: '{pdb_id}': no entity with chain id '{chain_id}' was found.\n\n") #todo log
    elif len(entities) > 1:
        #sys.stderr.write(f"Warning: '{pdb_id}': more than one entities with chain id '{chain_id}' were found:\n")
        #for x in entities:
        #        sys.stderr.write(f"{x[0]}\n")
        #sys.stderr.write(f"\n") #todo
        return entities
    else:
        return entities

#returns uniprotID, entityID, start, end
def get_uniprot_segments(pdb_id, chain_id): #todo cache
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}'

    try:
        parsedResponse = restAPI_get_json(url) #todo osetrit 404 error not found, napsat, ze nenalezen
    except HTTPError as err:
        if err.code == 404:
            raise ValueError(f"HTTP Error 404 - Not Found: structure {pdb_id} {chain_id} probably does not exist in UniProt or there is a problem with connection to the server.")
        else:
            raise
    uniprotRecords = parsedResponse[f"{pdb_id}"]["UniProt"]

    import re
    pattern = re.compile("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$") #UniProt accession format
    segments = []
    for id in uniprotRecords:
        uniprot_ID = id
        mappings = uniprotRecords[uniprot_ID]["mappings"]
        for entity in mappings:
            if entity["chain_id"] == chain_id:
                start_res_num = entity["start"]["residue_number"]
                end_res_num = entity["end"]["residue_number"]
                unp_start = entity["unp_start"]
                unp_end = entity["unp_end"]
                if (unp_end - unp_start != end_res_num - start_res_num):
                    raise ValueError(
                        f"Incorrect mapping: {pdb_id} {chain_id}, {uniprot_ID}: unp {unp_start}-{unp_end}, pdbe {start_res_num}-{end_res_num}")  # todo debug
                segments.append((uniprot_ID, unp_start, unp_end, start_res_num, end_res_num))

    if len(segments) < 1:
        raise ValueError(f"Structure {pdb_id} {chain_id} does not exist in UniProt. No uniprot segments were found. ")
    else:
        return segments

#returns parsed .json response (dictionary)
def restAPI_get(url):
    #todo check status?
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    return responseBody

#returns parsed .json response (dictionary)
def restAPI_get_json(url):
    return json.loads(restAPI_get(url))

#returns parsed .xml response (dictionary)
def restAPI_get_xml(url):
    return xmltodict.parse(restAPI_get(url)) #todo neni neco rychlejsiho? nepotrebuju ordered dictionary...

def get_fasta_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/FASTA/{pdb_id}{chain_id}.fasta"

def get_pdb_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/PDB/{pdb_id}{chain_id}.pdb"

def get_pdb_path2(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.pdb"

def get_lbs_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/lbs/{pdb_id}{chain_id}.txt"

def get_lbs_path2(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.txt"

def get_sasa_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/sasa/{pdb_id}{chain_id}.txt"

def get_mappings_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/mappings/{pdb_id}{chain_id}.txt"

def get_mappings_path2(data_dir, pdb_id, chain_id):
    return f"{data_dir}/{pdb_id}{chain_id}.txt"

def get_feature_path(data_dir, feature, pdb_id, chain_id):
    return f"{data_dir}/features/{feature}/{pdb_id}{chain_id}.txt"

def get_entity_id_old(pdb_id, chain_id):
    url = f"https://www.rcsb.org/pdb/rest/describeMol?structureId={pdb_id}.{chain_id}"
    response = str(restAPI_get(url))
    import re
    matches = re.findall("entityNr=\".+?\"", response)
    if (len(matches) < 1):
        raise ValueError(f"Structure {pdb_id} does not exist or does not have chain {chain_id}.")
    elif (len(matches) > 1):
        raise ValueError(f"Wrong number of matches: {matches}") #todo smazat, jen pro debugovani
    match=matches[0]
    entity_id = match.split('=')[1]
    entity_id = entity_id[1:-1] #remove ""
    return int(entity_id)

def get_entity_id(pdb_id, chain_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
    response = restAPI_get_json(url)
    molecules = response[pdb_id]
    entity_ids = []
    for molecule in molecules:
        if "polypeptide" in molecule["molecule_type"]:
            if chain_id in molecule["in_chains"]:
                entity_ids.append(molecule["entity_id"])
    if len(entity_ids) != 1:
        raise ValueError(f"{pdb_id} {chain_id}: Wrong number of entity IDs! {len(entity_ids)}") #todo jen debug
    #todo kdyz neexistuje ten chain, vrati to 0 entity IDS -> nejaka hlaska
    return int(entity_ids[0]) #todo

def res_mappings_author_to_pdbe(pdb_id, chain_id, cache_file=""):
    if (cache_file != ""):
        mappings = np.genfromtxt(cache_file, delimiter=' ', dtype='str')
        return list(mappings)
    else:
        mappings = []
        response = restAPI_get_json(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id}/chain/{chain_id}")
        entity_id = get_entity_id(pdb_id, chain_id)
        molecules = response[pdb_id]["molecules"]
        count = 0
        for molecule in molecules:
            if molecule["entity_id"] == entity_id:
                for residue in molecule["chains"][0]["residues"]:
                    key = str(residue["author_residue_number"]) + residue["author_insertion_code"] #author residue number
                    val = residue["residue_number"] #pdbe molecule residue number
                    mappings.append((key, val))
                count += 1
        if count != 1:
            raise ValueError(f"Error: More or less than one molecule with entity number {entity_id} was found.")
            return
        return list(mappings)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            line = line.split()
            pdb_id = line[0]
            chain_ids = line[1] #todo check
            chain_ids = chain_ids.upper()
            list.append((pdb_id, chain_ids))
    return list

def parse_dataset_split_chains(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            line = line.split()
            pdb_id = line[0]
            chain_ids = line[1]  # todo check
            chain_ids = chain_ids.upper()
            for chain in chain_ids.split(','):
                list.append((pdb_id, chain))
    return list

def parse_prank_dataset(filepath):
    pdb_files = []
    with open(filepath) as f:
        for line in f:
            #todo jaka je presna struktura .ds souboru v pranku? nechybi mi neco?
            if (line[0] != '#'): #skip comments
                words = line.split()
                if (len(words) > 0 and words[0] != 'HEADER:'): # skip empty lines and header
                    pdb_files.append(words[0])

    return pdb_files

def getCenterOfMass(atoms):
    totalMass = 0.0
    x = 0
    y = 0
    z = 0
    for a in atoms:
        m = a.mass
        totalMass += m
        x += a.coord[0] * m
        y += a.coord[1] * m
        z += a.coord[2] * m
    return (x / totalMass, y / totalMass, z / totalMass)

def getFullAuthorResNum(residue_id):
    auth_res_num = str(residue_id[1])
    if not (residue_id[2].isspace()): # biopython returns a space instead of empty string
        auth_res_num += str(residue_id[2])  # insertion code
    return auth_res_num

def isPartOfChain(residue, mappings):
    if (residue.id[0].isspace()):  # hetero flag is empty
        return True
    elif (is_aa(residue, standard=False)):  # HETATM, but nonstandard AA code (MSE, LYZ etc.)
        # check if mapping exist for this residue number:
        auth_res_num = getFullAuthorResNum(residue.id)
        if (auth_res_num in mappings):
            return True # is part of the chain, not ligand
        else:
            return False  # is truly a ligand
    else:
        return False

def getStructuresFromDirectory(dir):
    list = []
    for filename in os.listdir(dir):
        pdb_id = filename[:4]
        chain_id = filename[4:5]
        list.append((pdb_id, chain_id))

    return list
