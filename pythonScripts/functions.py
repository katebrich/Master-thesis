import urllib.request
import json
import sys

#returns uniprotID, entityID, start, end
def get_uniprot_entity(pdb_id, chain_id):
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'
    parsedResponse = restAPI_get(url)
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
    return json.loads(responseBody)

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

def get_mmcif_path(data_dir, pdb_id, chain_id):
    return f"{data_dir}/mmCIF/{pdb_id}{chain_id}.cif"
