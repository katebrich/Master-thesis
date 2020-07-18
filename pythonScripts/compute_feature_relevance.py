from get_ligand_binding_sites import get_ligand_binding_sites
from get_feature import get_feature
from parse_dataset import parse_dataset
from functions import get_pdb_path
import statistical_analysis as stats
import pandas as pd
import os
import pickle

dataset_name = "coach420"
dataset_path = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'
data_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}/" #todo parametr

dataset = parse_dataset(dataset_path)
feature = "dynamine" #todo
pairs = []
i = 1
total = len(dataset)

for structure in dataset:
    #todo
    pdb_id = structure[0]
    chain_id = structure[1]

    print(f"Processing structure {pdb_id} {chain_id}")

    #todo definvat co ocekavam za typy

    cache_file = f"{data_dir}lbs/{pdb_id}{chain_id}.txt"
    if os.path.isfile(cache_file): # read cached values
        with open(cache_file, "rb") as fp:
            lbs = pickle.load(fp)
    else:
        lbs = get_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_id))
        with open(cache_file, "wb") as fp:  # save cache
            pickle.dump(lbs, fp)

    feature_vals = get_feature(feature, data_dir, pdb_id, chain_id)

    #pair feature values with ligand binding sites
    # numbering is corresponding to the PDBe molecule
    missing_vals = []
    lbs_dict = dict(lbs)
    for val in feature_vals:
        res_num = val[0]
        feature_val = val[1]
        if res_num in lbs_dict:
            lbs_val = lbs_dict[res_num]
        else:
            missing_vals.append(res_num)
            continue;
        pairs.append((lbs_val, feature_val))

    if (len(missing_vals) > 0):
        print(f"Warning: missing feature values for residues: {missing_vals}")
    print(f"{i}/{total}: structure {pdb_id} {chain_id} processed.")
    i += 1


#take the whole dataset and run statistical analysis
#print(pairs)

#stats.compute_AA_frequencies(pairs)

stats.fischers_exact_test(pairs)
