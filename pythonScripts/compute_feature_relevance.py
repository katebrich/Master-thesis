from get_ligand_binding_sites import get_ligand_binding_sites
from get_feature import get_PTM
import statistical_analysis as stats
import pandas as pd

dataset = [("2rbt","X"), ("3cx5","O")] #todo
feature = "PTM"

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    lbs = get_ligand_binding_sites(pdb_id, chain_id, f'/home/katebrich/Documents/diplomka/PDBe_files/pdb{structure[0]}.ent')
    feature_vals = get_PTM(pdb_id, chain_id)
    #feature_vals[0] = 1 #todo smazat
    #todo sanity check - same length
    if (len(lbs) != len(feature_vals)):
        print("error - length of lbs and feature_vals vectors is not the same")

    df = pd.DataFrame(list(zip(lbs, feature_vals)),
                      columns=['ligand_binding_sites', "feature"])

    stats.fischers_exact_test(df)
