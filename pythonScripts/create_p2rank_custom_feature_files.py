import sys
import getopt
from Bio.PDB import PDBParser

from Config import Config
from helper import parse_dataset, get_pdb_path
import os
import pandas as pd
from helper import res_mappings_author_to_pdbe
import csv
import numpy as np
import Logger
import time
from multiprocessing import Pool, Value

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class P2RankCustomFeatureCreator():
    dataset_file = ""
    output_dir = ""
    input_dir = ""
    features = []
    defaults = []
    total = ""
    mappings_cache_dir=""

    def __init__(self, dataset_file, output_dir, input_dir, features):
        self.dataset_file = dataset_file
        self.output_dir = output_dir
        self.input_dir = input_dir
        self.features = features
        self.mappings_cache_dir = os.path.join(input_dir, "mappings")  # todo

        config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.json")
        config = Config(config_path)
        for feature in features:
            self.defaults.append(config.get_feature_default(feature)) #todo test

    def run(self, threads=1):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset(dataset_file)

        start = time.time()
        logger.info(f"Creating p2rank custom feature files started...Features: {self.features}") #todo log asi rovnou do konzole, kdyz se to spousti takhle samostatne

        self.total = len(dataset)

        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.create_file, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Creating p2rank custom feature files finished: All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Creating p2rank custom feature files finished: Some structures were not processed successfully: \n{errors_format}")
        logger.debug(f"Finished in {time.time() - start}")

    def read_feature_vals(self, feature, data_dir, pdb_id, chain_id):
        path = os.path.join(data_dir, "features", feature, f"{pdb_id}{chain_id}.txt")  # todo
        feature_vals = np.genfromtxt(path, delimiter=' ', dtype='str')
        return dict(feature_vals)

    def create_file(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        errors = []
        try:
            aa_codes = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]  # TODO smazat az se odstrani chyba v pranku

            df = pd.DataFrame(columns=["chain", "ins. code", "seq. code"] + features)
            parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
            pdb_path = get_pdb_path(input_dir, pdb_id, chain_id)
            struct = parser.get_structure(pdb_id, pdb_path)

            feature_vals = []
            for feature in features:
                feature_vals.append(self.read_feature_vals(feature, input_dir, pdb_id, chain_id))  # todo

            mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, os.path.join(self.mappings_cache_dir,
                                                                                       f"{pdb_id}{chain_id}.txt")))  # todo helper
            for residue in struct.get_residues():
                if (residue.id[0][2:] in aa_codes):  # todo smazat fix pro prank
                    if (residue.id[2].isspace()):
                        ins_code = ""
                    else:
                        ins_code = residue.id[2]
                    seq_code = residue.id[1]
                    # res_num = str(seq_code) + str(ins_code)
                    feat_tuple = tuple(self.defaults)
                    df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code,
                                                                                      seq_code) + feat_tuple
                elif (residue.id[0].isspace() or residue.id[
                    0] == "H_MSE"):  # skip hetero-residues #todo selenomethionine??? #todo nonstandard residues...jako pri pocitani LBS
                    if (residue.id[2].isspace()):  # biopython returns a space instead of empty string
                        ins_code = ""
                    else:
                        ins_code = residue.id[2]
                    seq_code = residue.id[1]
                    auth_res_num = str(seq_code) + str(ins_code)
                    pdbe_res_num = mappings[auth_res_num]
                    feat_list = []
                    # todo missing feature vals for res_num
                    for j in range(0, len(features)):
                        val = feature_vals[j].get(pdbe_res_num, self.defaults[j])
                        feat_list.append(val)
                    feat_tuple = tuple(feat_list)
                    df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code,
                                                                                      seq_code) + feat_tuple
            output_path = os.path.join(output_dir, os.path.basename(pdb_path) + ".csv")
            with open(output_path, 'w') as file:
                file.write(df.to_csv(index=False, quoting=csv.QUOTE_ALL))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append(structure)
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED !")
            else:
                logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    # todo zrychlit to!!!

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args





dataset_file = ""
output_dir = ""
input_dir = ""
features = ""
threads = 1

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:i:t:')
except getopt.GetoptError as err:
    logger.error(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-i", "--input_dir"): #folder with PDB files
        input_dir = arg
    elif opt in ("-t", "--threads"):
        threads = arg #todo check if threads >= 1, int
    elif opt in ("-f", "--features"):
        features = arg.split(',')

if (dataset_file == ""):
    logger.error("Dataset must be specified.")
    sys.exit(1)
if (output_dir == ""):
    logger.error("Output directory must be specified.")
    sys.exit(1)
if (input_dir == ""):
    logger.error("Input directory must be specified.")
    sys.exit(1)
#todo check features

p = P2RankCustomFeatureCreator(dataset_file, output_dir, input_dir, features)
p.run(threads)


