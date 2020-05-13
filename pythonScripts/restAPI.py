import requests, sys

accessions = {"Q14676", "P69892"}

#modified residues for protein Q14676
requestURL = f"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession={','.join(accessions)}&types=MOD_RES"
r = requests.get(requestURL, headers={ "Accept" : "application/json"})
if not r.ok:
  r.raise_for_status()
  sys.exit()
responseBody = r.text
print(responseBody)

import json

#returns dictionary
parsedResponse = json.loads(responseBody)

#iterate over all objects (=Uniprot accessions)
for i in range(0, len(accessions)):
    features = parsedResponse[i]["features"]
    #get list of modified residues
    mod_res = []
    for feature in features:
        mod_res.append(feature["begin"])
    print(mod_res)
