import urllib.request
import json
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'p:c:', ["pdb_id", "chain_id"])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    #todo print help
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        # todo print help
        sys.exit()
    elif opt in ("-p", "--pdb_id"):
        pdb_id = arg
    elif opt in ("-c", "--chain_id"):
        chain_id = arg
    #todo else unknown option

url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'

req = urllib.request.Request(url)
with urllib.request.urlopen(req) as f:
 responseBody = f.read()


#returns dictionary of dictionaries
parsedResponse = json.loads(responseBody)

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
    print(f"Error: '{pdb_id}': no entity with chain id '{chain_id}' was found.")
    sys.exit(1)
elif len(entities) > 1:
    print(f"Error: '{pdb_id}': more than one entities with chain id '{chain_id}' were found:")
    for x in entities:
            print(x[0])
    sys.exit(1)
else:
    print(entities[0])
    sys.exit(0)