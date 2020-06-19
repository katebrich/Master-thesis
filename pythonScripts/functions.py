import urllib.request
import json
import sys

#returns uniprotID, entityID, start, end
def get_uniprot_entity(pdb_id, chain_id):
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'
    parsedResponse = restAPI_get(url)
    uniprotRecords = parsedResponse[f"{pdb_id}"]["UniProt"]

    #find entity with given chain_id
    entities = []
    for record in uniprotRecords:
        uniprot_ID = record
        mappings = uniprotRecords[record]["mappings"]
        for entity in mappings:
            if entity["chain_id"] == chain_id:
                entities.append((uniprot_ID, entity["entity_id"], entity["unp_start"], entity["unp_end"]))

    if len(entities) < 1:
        sys.stderr.write(f"Error: '{pdb_id}': no entity with chain id '{chain_id}' was found.\n\n")
    elif len(entities) > 1:
        sys.stderr.write(f"Error: '{pdb_id}': more than one entities with chain id '{chain_id}' were found:\n")
        for x in entities:
                sys.stderr.write(f"{x[0]}\n")
        sys.stderr.write(f"\n")
    else:
        return entities[0]

#returns parsed .json response (dictionary)
def restAPI_get(url):
    #todo check status
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    return json.loads(responseBody)