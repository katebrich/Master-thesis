import getopt
import os
import sys
import numpy as np

from helper import parse_dataset, eprint
from statistical_analysis import welchs_t_test, fischers_exact_test
from features import types_of_features

dataset_file = ""
output_dir = ""
lbs_dir = ""
feature_dir = ""
feature_name = ""
threads = 1 #todo asi nebude potreba

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:v:d:o:l:t:')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-l", "--lbs_dir"):
        lbs_dir = arg
    elif opt in ("-t", "--threads"):
        threads = arg #todo check if threads >= 1, int
    elif opt in ("-v", "--feature_values_dir"):
        feature_dir = arg
    elif opt in ("-f", "--feature_name"):
        feature_name = arg

if (dataset_file == ""):
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (dataset_file == ""):
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (feature_name == ""):
    eprint("ERROR: Feature name must be specified.")
    sys.exit(1)
if (lbs_dir == ""):
    eprint("ERROR: Directory with ligand binding sites data must be specified.")
    sys.exit(1)
if (feature_dir == ""):
    eprint("ERROR: Directory with feature values must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#todo else smazat??


dataset = parse_dataset(dataset_file)
pairs = []
total = len(dataset)
i = 1
#pair feature values with ligand binding sites
for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]
    # numbering is corresponding to the PDBe molecule
    missing_vals = []

    # get ligand binding sites values
    output_file = f"{lbs_dir}/{pdb_id}{chain_id}.txt"
    lbs = np.genfromtxt(output_file, delimiter=' ')
    lbs_dict = dict(lbs)

    # get feature values
    output_file = f"{feature_dir}/{pdb_id}{chain_id}.txt"
    feature = np.genfromtxt(output_file, delimiter=' ')
    feature_vals = list(feature)

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

type_of_feature = types_of_features[feature_name]

if (type_of_feature == "discrete"):
    fischers_exact_test(pairs)
elif (type_of_feature == "continuous"):
    welchs_t_test(pairs)
else:
    eprint(f"ERROR: Unknown type of feature '{type_of_feature}'")
    sys.exit(1)
