from helper import *

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

class SecondaryStructure():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/secondary_structure/{pdb_id}/{entity_id}"
        response = restAPI_get_json(url)

        feature_vector = ['X'] * int(response[pdb_id]["length"])
        for rec in response[pdb_id]["data"]:
            dataType = rec["dataType"]  # Helix, Strand, MobiDB
            if (dataType == "Helix"):
                val = 'H'
            elif (dataType == "Strand"):
                val = 'S'
            elif (dataType == "MobiDB"):  # skip this type
                continue
            else:
                raise ValueError(f"PDBe KB - secondary structure: {pdb_id} {chain_id} : unknown dataType {dataType}")
            for resi in rec["residues"]:
                feat_begin = int(resi["startIndex"])
                feat_end = int(resi["endIndex"])
                for i in range(feat_begin, feat_end + 1):
                    feature_vector[i - 1] = val

        feature_vals = []
        for res_num in range(1, int(response[pdb_id]["length"] + 1)):
            feature_vals.append((res_num, feature_vector[res_num - 1]))

        return feature_vals

class Conservation():
    def get_values(self, data_dir, pdb_id, chain_id):
        entity_id = get_entity_id(pdb_id, chain_id)
        url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/sequence_conservation/{pdb_id}/{entity_id}"
        response = restAPI_get_json(url)

        feature_vals = []
        for resi in response[pdb_id]["data"]:
            feat_begin = int(resi["start"])
            feat_end = int(resi["end"])
            if (feat_begin != feat_end):
                raise ValueError("Pdb KB conservation: start is not same as end!!")  # todo only for debug
            feature_vals.append((feat_begin, resi["conservation_score"]))

        return feature_vals

class Dynamine():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("dynamine", "backbone", pdb_id, chain_id)

class Efoldmine():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("dynamine", "efoldmine", pdb_id, chain_id)

class Depth():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_funPDBe("depth", "complex_residue_depth", pdb_id, chain_id)