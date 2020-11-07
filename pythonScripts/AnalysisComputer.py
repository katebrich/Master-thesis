import os
import sys
import time

import numpy as np

from helper import getStructuresFromDirectory
from unused.statistical_analysis import welchs_t_test, fischers_exact_test, chi_squared_test
#from unused.features import types_of_features
import logger

logger = logger.get_logger(os.path.basename(__file__))

class AnalysisComputer():
    output_dir = ""
    lbs_dir = ""
    config = ""
    feature_dir = ""
    feature_name = ""
    errors = []
    pairs = []
    total = ""
    counter = 1

    def __init__(self, output_dir, lbs_dir, config):
        self.output_dir = output_dir
        self.lbs_dir = lbs_dir
        self.config = config

    def run(self, feature_name, feature_dir):
        self.feature_dir = feature_dir
        self.feature_name = feature_name

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = getStructuresFromDirectory(self.lbs_dir)

        start = time.time()
        logger.info(f"Running analysis of feature {self.feature_name} started...")

        self.pairs = []
        self.errors = []
        self.total = len(dataset)
        self.counter = 1

        #pair feature values with ligand binding sites
        for structure in dataset:
            self.compute_pairs(structure)

        if (len(self.errors) == 0):
            logger.info(f"All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in self.errors)
            logger.warning(f"Some structures were not processed successfully: skipping them in the analysis...\n{errors_format}")

        feature_type = self.config.get_feature_type(self.feature_name) #todo kdyz neznamy typ, tak aby to cele nespadlo..asi

        file = os.path.join(self.output_dir, f"{feature_name}.txt")
        if (feature_type == "binary"):
            logger.info("Running Fischer's exact test")
            fischers_exact_test(self.pairs, file)
        elif (feature_type == "continuous"):
            logger.info("Running Welch's T-test")
            welchs_t_test(self.pairs, file)
        elif (feature_type == "categorical"):
            logger.info("Running Chi-squared test")
            chi_squared_test(self.pairs, file)
        else:
            #todo
            logger.error(f"Unknown type of feature '{feature_type}'. Please specify the type in config.")
            sys.exit(1)

        logger.info(f"Running analysis finished. Results saved to {file}")

        logger.debug(f"Finished in {time.time() - start}")

    def compute_pairs(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        missing_vals = []

        try:
            # get ligand binding sites values
            file = os.path.join(self.lbs_dir, f"{pdb_id}{chain_id}.txt")
            lbs = np.genfromtxt(file, delimiter=' ', dtype=None)
            lbs_dict = dict(lbs)

            # get feature values
            file = os.path.join(self.feature_dir, f"{pdb_id}{chain_id}.txt")
            feature = np.genfromtxt(file, delimiter=' ', dtype=None, encoding=None)
            feature_vals = list(feature)

            for val in feature_vals:
                res_num = val[0]
                feature_val = val[1]
                if res_num in lbs_dict:
                    lbs_val = lbs_dict[res_num]
                else:
                    missing_vals.append(res_num)
                    continue
                self.pairs.append((lbs_val, feature_val))
            if (len(missing_vals) > 0):
                logger.debug(f"{pdb_id} {chain_id}: Missing feature values for residues: {missing_vals}") #todo vyresit to
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            self.counter  += 1
            if (error):
                self.errors.append(structure)
                logger.error(f"{self.counter}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED !")
            else:
                logger.debug(f"{self.counter}/{self.total}: {pdb_id} {chain_id} processed")


'''
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
'''
