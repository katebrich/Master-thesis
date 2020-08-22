import getopt
import os
import sys
import traceback

import numpy as np

from helper import parse_dataset_split_chains
from statistical_analysis import welchs_t_test, fischers_exact_test, chi_squared_test
from features import types_of_features
import logger

logger = logger.get_logger(os.path.basename(__file__))

def compute_pairs(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error = False
    missing_vals = []

    try:
        # get ligand binding sites values
        file = os.path.join(lbs_dir, f"{pdb_id}{chain_id}.txt")
        lbs = np.genfromtxt(file, delimiter=' ', dtype=None)
        lbs_dict = dict(lbs)

        # get feature values
        file = os.path.join(feature_dir, f"{pdb_id}{chain_id}.txt")
        feature = np.genfromtxt(file, delimiter=' ', dtype=None)
        feature_vals = list(feature)

        for val in feature_vals:
            res_num = val[0]
            feature_val = val[1]
            if res_num in lbs_dict:
                lbs_val = lbs_dict[res_num]
            else:
                missing_vals.append(res_num)
                continue;
            global pairs
            pairs.append((lbs_val, feature_val))
        if (len(missing_vals) > 0):
            logger.debug(f"Missing feature values for residues: {missing_vals}") #todo vyresit to
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as ex:
        error = True
        logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}")
    finally:
        global counter
        counter  += 1
        if (error):
            errors.append(structure)
            logger.error(f"{counter}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
        else:
            logger.debug(f"{counter}/{total}: {pdb_id} {chain_id} processed")

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
    logger.error(err) #unknown option or missing argument
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
    logger.error("Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (dataset_file == ""):
    logger.error("Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (feature_name == ""):
    logger.error("Feature name must be specified.")
    sys.exit(1)
if (lbs_dir == ""):
    logger.error("Directory with ligand binding sites data must be specified.")
    sys.exit(1)
if (feature_dir == ""):
    logger.error("Directory with feature values must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#todo else smazat??

dataset = parse_dataset_split_chains(dataset_file)

logger.info("Running analysis started...")

pairs = []
total = len(dataset)
counter = 1
errors = []
#pair feature values with ligand binding sites
for structure in dataset:
    compute_pairs(structure)

if (len(errors) == 0):
    logger.info(f"All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in errors)
    logger.warning(f"Some structures were not processed successfully: skipping them in the analysis...\n{errors_format}")


#stats.compute_AA_frequencies(pairs)

type_of_feature = types_of_features[feature_name]

file = os.path.join(output_dir, "results.txt")
if (type_of_feature == "binary"):
    logger.info("Running Fischer's exact test")
    fischers_exact_test(pairs, file)
elif (type_of_feature == "continuous"):
    logger.info("Running Welch's T-test")
    welchs_t_test(pairs, file)
elif (type_of_feature == "categorical"):
    logger.info("Running Chi-squared test")
    chi_squared_test(pairs, file)
else:
    logger.error(f"Unknown type of feature '{type_of_feature}'")
    sys.exit(1)

logger.info(f"Running analysis finished. Results saved to {file}")