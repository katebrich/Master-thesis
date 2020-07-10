from get_ligand_binding_sites import get_ligand_binding_sites
from get_feature import get_feature
from parse_dataset import parse_dataset
from functions import get_mmcif_path
import statistical_analysis as stats
import pandas as pd

dataset_path = f'/home/katebrich/Documents/diplomka/datasets/test.ds'
data_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420" #todo podle datasetu

dataset = parse_dataset(dataset_path)
feature = "hydropathy" #todo

df = pd.DataFrame(columns = ["ligand_binding_sites", "feature"])

for structure in dataset:
    #todo
    pdb_id = structure[0]
    chain_id = structure[1]

    print(f"Processing structure {pdb_id} {chain_id}")

    #todo definvat co ocekavam za typy
    lbs = get_ligand_binding_sites(pdb_id, chain_id, get_mmcif_path(data_dir, pdb_id, chain_id)) #todo
    feature_vals = get_feature(feature, data_dir, pdb_id, chain_id)

    #pair feature values with ligand binding sites
    # numbering is corresponding to the PDBe molecule
    pairs = []
    missing_vals = []
    lbs_dict = dict(lbs)
    for val in feature_vals:
        res_num = val[0]
        feature_val = val[1]
        if str(res_num) in lbs_dict: #todo proc v dictionary string, ale ve feature_vals int? opravit
            lbs_val = lbs_dict[str(res_num)]
        else:
            missing_vals.append(res_num)
            continue;
        pairs.append((lbs_val, feature_val))

    new_df = pd.DataFrame(pairs,
                     columns=['ligand_binding_sites', "feature"])

    df = pd.concat([df, new_df])

    if (len(missing_vals) > 0):
        print(f"Warning: missing feature values for residues: {missing_vals}")
    print(f"structure {pdb_id} {chain_id} processed.")
#take the whole dataset and run statistical analysis
stats.fischers_exact_test(df)
