from helper import get_uniprot_entity, restAPI_get_json
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'p:')
except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    #todo print help
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        # todo print help
        sys.exit()
    elif opt in ("-p"):
        pdb_id = arg[:4]
        chain_id = arg[4]
    #todo else unknown option

#pdb_id = "1a2k"
#chain_id = "C"

try:
    uniprot_id, entity_id, begin, end = get_uniprot_entity(pdb_id, chain_id)
except:
    sys.exit(3)

print(pdb_id)
print(chain_id)
print(uniprot_id)
print(entity_id)
print(begin)
print(end)

# pdbe-kb - vraci neco uplne jineho nez uniprot a navic ani ne to same jako je na webu pdbe-kb...
# url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/modified_AA_or_NA/{pdb_id}"

#uniprot - jen MOD_RES type (podkategorie PTM)
#url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types=MOD_RES"

url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=PTM"

response = restAPI_get_json(url)

feature_vector = [0] * (end-begin+1) #including both start and end AAs

for feature in response["features"]:
    feat_begin = int(feature["begin"])
    feat_end = int(feature["end"])
    type = feature["type"]
    if (type == "DISULFID"): #disulfide bond - 'begin' and 'end' mean connected AAs
        rng = {feat_begin, feat_end}
    else:
        rng = range(feat_begin, feat_end+1)
    for i in rng:
        res = i - begin # mapping pdb residues to uniprot entry, counting from 0
        if res >= 0 and res < end-begin+1:
            feature_vector[res] = 1

print(*feature_vector, sep=",")
#print(response)
print()



