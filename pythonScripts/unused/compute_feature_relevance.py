from compute_ligand_binding_sites import __compute_ligand_binding_sites
from features import get_feature
from helper import parse_dataset
from helper import get_pdb_path
from statistical_analysis import welchs_t_test, fischers_exact_test
import pandas as pd
import os
import pickle

cache = True
dataset_name = "chen11"
dataset_path = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'
data_dir = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/" #todo parametr

dataset = parse_dataset(dataset_path)
feature = "unp_PTM" #todo
pairs = []
i = 1
total = len(dataset)

if (cache):
    lbs_cache_dir = f"{data_dir}lbs/"
    if not os.path.exists(lbs_cache_dir):
        os.makedirs(lbs_cache_dir)

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    print(f"Processing structure {pdb_id} {chain_id}")

    if (cache):
        cache_file = f"{lbs_cache_dir}/{pdb_id}{chain_id}.txt"
        if os.path.isfile(cache_file): # read cached values
            with open(cache_file, "rb") as fp:
                lbs = pickle.load(fp)
        else:
            lbs = __compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_id))
            with open(cache_file, "wb") as fp:  # save cache
                pickle.dump(lbs, fp)
    else:
        lbs = __compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_id))

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

#stats.fischers_exact_test(pairs)
stats.welchs_t_test(pairs)