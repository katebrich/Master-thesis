from functions import get_uniprot_entity, restAPI_get_json

pdb_id = "3cx5"
chain_id = "O"

uniprot_id, entity_id, start, end = get_uniprot_entity(pdb_id, chain_id)

url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/secondary_structure/{pdb_id}/{entity_id}"

response = restAPI_get_json(url)

print(response)



