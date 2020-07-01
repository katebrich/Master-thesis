from functions import get_uniprot_entity, restAPI_get

def get_PTM(pdb_id, chain_id):
    try:
        entities = get_uniprot_entity(pdb_id, chain_id)
    except:
        print("error")
        #sys.exit(3)
        #todo

    feature_vals = []

    for entity in entities:
        uniprot_id = entity[0]
        entity_id = entity[1]
        unp_start = entity[2]
        unp_end = entity[3]
        start_res_num = entity[4]
        end_res_num = entity[5]

        #print(pdb_id)
        #print(chain_id)
        #print(uniprot_id)
        #print(entity_id)
        #print(unp_start)
        #print(unp_end)
        #print(start_res_num)
        #print(end_res_num)

        # pdbe-kb - vraci neco uplne jineho nez uniprot a navic ani ne to same jako je na webu pdbe-kb...
        # url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/modified_AA_or_NA/{pdb_id}"

        # uniprot - jen MOD_RES type (podkategorie PTM)
        # url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?types=MOD_RES"

        url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}?categories=PTM"

        response = restAPI_get(url)

        feature_vector = [0] * (unp_end - unp_start + 1)  # including both start and end AAs

        for feature in response["features"]:
            feat_begin = int(feature["begin"])
            feat_end = int(feature["end"])
            type = feature["type"]
            if (type == "DISULFID"):  # disulfide bond - 'begin' and 'end' mean connected AAs
                rng = {feat_begin, feat_end}
            else:
                rng = range(feat_begin, feat_end + 1)
            for i in rng:
                res = i - unp_start  # mapping pdb residues to uniprot entry, counting from 0
                if res >= 0 and res < unp_end - unp_start + 1:
                    feature_vector[res] = 1

        i = 0
        for res_num in range(start_res_num, end_res_num + 1):
            feature_vals.append((res_num, feature_vector[i]))
            i += 1

        #print(*feature_vector, sep=",")
        #print(response)

    return feature_vals