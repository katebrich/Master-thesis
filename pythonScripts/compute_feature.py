import getopt
import itertools
import os
import sys
import threading
import traceback
import time

from helper import eprint, parse_dataset_split_chains
from features import *
import logger

logger = logger.get_logger(os.path.basename(__file__))

def __compute_feature(structure, name_of_feature):
    pdb_id = structure[0]
    chain_id = structure[1]
    error=False
    try:
        feat_vals = get_feature(name_of_feature, input_dir, pdb_id, chain_id)

        output_file = f"{output_dir}/{pdb_id}{chain_id}.txt"
        with open(output_file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in feat_vals))
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as ex:
        error=True
        logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}", exc_info=True)
    finally:
        with threadLock:
            global counter
            idx = counter
            counter += 1
        if (error):
            errors.append(structure)
            logger.error(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
        else:
            logger.debug(f"{idx}/{total}: {pdb_id} {chain_id} processed")


dataset_file = ""
output_dir = ""
input_dir = ""
threads = 1
feature=""

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
    elif opt in ("-f", "--feature"):
        feature = arg

if (dataset_file == ""):
    logger.error("Dataset must be specified.")
    sys.exit(1)
if (output_dir == ""):
    logger.error("Output directory must be specified.")
    sys.exit(1)
if (input_dir == ""):
    logger.error("Input directory must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#todo else smazat??

dataset = parse_dataset_split_chains(dataset_file) #todo co kdyz neni spravny format

start = time.time()
logger.info(f"Computing feature {feature} started...")

total = len(dataset)
counter = 1
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur
errors=[]

# run dynamine in 1 thread only !!!
if (feature == "dynamine"): #todo neni asi potreba
    threads = 1

if (threads == 1):
    for structure in dataset:
        __compute_feature(structure, feature)
else:
    args = zip(dataset, itertools.repeat(feature, len(dataset)))
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.starmap(__compute_feature, args)
    pool.close()

if (len(errors) == 0):
    logger.info(f"Computing feature {feature} finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in errors)
    logger.warning(f"Computing feature {feature} finished: Some structures were not processed successfully: \n{errors_format}")
logger.debug(f"Finished in {time.time() - start}")