from __future__ import print_function
import urllib.request
import json
import sys
import xmltodict


#returns uniprotID, entityID, start, end
def get_uniprot_entity(pdb_id, chain_id):
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'
    parsedResponse = restAPI_get_json(url)
    uniprotRecords = parsedResponse[f"{pdb_id}"]["UniProt"]

    #find entity with given chain_id
    entities = [] #todo dictionary?
    for record in uniprotRecords:
        uniprot_ID = record
        mappings = uniprotRecords[record]["mappings"]
        for entity in mappings:
            start_res_num = entity["start"]["residue_number"]
            end_res_num = entity["end"]["residue_number"]
            if entity["chain_id"] == chain_id:
                entities.append((uniprot_ID, entity["entity_id"], entity["unp_start"], entity["unp_end"], start_res_num, end_res_num))

    if len(entities) < 1:
        sys.stderr.write(f"Error: '{pdb_id}': no entity with chain id '{chain_id}' was found.\n\n")
    elif len(entities) > 1:
        sys.stderr.write(f"Warning: '{pdb_id}': more than one entities with chain id '{chain_id}' were found:\n")
        for x in entities:
                sys.stderr.write(f"{x[0]}\n")
        sys.stderr.write(f"\n")
        return entities
    else:
        return entities

#returns parsed .json response (dictionary)
def restAPI_get(url):
    #todo check status
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

def isInDistance(threshold, residue1, residue2):
    #todo optimize
    for atom1 in residue1.child_list:
        for atom2 in residue2.child_list:
            if atom1 - atom2 < threshold: # - operator measures distance
                #print(residue1.id[1])
                return True
    return False

def get_fasta_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/FASTA/{pdb_id}{chain_id}.fasta"

def get_pdb_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/PDB/{pdb_id}{chain_id}.pdb"

#def get_mmcif_path(data_dir, pdb_id, chain_id):
#    return f"{data_dir}/mmCIF/{pdb_id}{chain_id}.cif"

def get_entity_id(pdb_id, chain_id):
    url = f"https://www.rcsb.org/pdb/rest/getEntityInfo?structureId={pdb_id}"
    parsedResponse = restAPI_get_xml(url)
    entities = parsedResponse["entityInfo"]["PDB"]["Entity"]
    entities_list = []
    if (type(entities) is not list): #jen 1 entita, rovnou ty atributy
        entities_list.append(entities)
    else:
        entities_list = entities
    if len(entities_list) == 0:
        print(f"Error: no entity found for {pdb_id} {chain_id}") #todo tohle asi nebude fungovat
        return
    entity_id = ""
    for entity in entities_list:
        success=False
        chains_list = []
        chains = entity["Chain"]
        if (type(chains) is not list):  # jen 1 entita, rovnou ty atributy
            chains_list.append(chains)
        else:
            chains_list = chains
        for chain in chains_list:
            if (chain["@id"] == chain_id):
                success=True
                break
        if success:
            new_entity_id = entity["@id"]
            if new_entity_id != entity_id and entity_id != "":
                print(f"Error: not all the entity IDs for {pdb_id} {chain_id} are the same!")
                return #todo tohle pak smazat
            entity_id = new_entity_id
    return int(entity_id)

def res_mappings_author_to_pdbe(pdb_id, chain_id): #todo cache
    mappings = []
    response = restAPI_get_json(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id}/chain/{chain_id}")
    entity_id = get_entity_id(pdb_id, chain_id)
    molecules = response[pdb_id]["molecules"]
    count = 0
    for molecule in molecules:
        if molecule["entity_id"] == entity_id:
            for residue in molecule["chains"][0]["residues"]:
                key = str(residue["author_residue_number"]) + residue["author_insertion_code"] #author residue number
                val = residue["residue_number"] #pdbe residue number
                mappings.append((key, val))
            count += 1
    if count != 1:
        print(f"Error: More or less than one molecule with entity number {entity_id} was found.")
        return
    return dict(mappings)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            pdb_id, chain_ids = line.split()
            chain_ids = chain_ids.upper()
            list.append((pdb_id, chain_ids))
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