from get_ligand_binding_sites import get_ligand_binding_sites
from get_feature import get_PTM
from parse_dataset import parse_dataset
import statistical_analysis as stats
import pandas as pd

filepath = f'/home/katebrich/Documents/diplomka/datasets/test.ds'

dataset = parse_dataset(filepath)
feature = "PTM" #todo

df = pd.DataFrame(columns = ["ligand_binding_sites", "feature"])

for structure in dataset:
    #todo
    pdb_id = structure[0]
    chain_id = structure[1]

    lbs = get_ligand_binding_sites(pdb_id, chain_id, f'/home/katebrich/Documents/diplomka/PDBe_files/coach420/{pdb_id}{chain_id}.pdb')
    feature_vals = get_PTM(pdb_id, chain_id)

    #pair feature values with ligand binding sites
    # !! feature_vals numbered from 1 -> for example feature_vals[2] pairs with lbs[1]
    pairs = []
    for val in feature_vals:
        pairs.append((lbs[val[0] - 1], val[1]))

    new_df = pd.DataFrame(pairs,
                     columns=['ligand_binding_sites', "feature"])

    df = pd.concat([df, new_df])

    print(f"structure {pdb_id} {chain_id} processed.")

print(df)
#take the whole dataset and run statistical analysis
stats.fischers_exact_test(df)
