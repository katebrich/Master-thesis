import getopt
import itertools
import os
import sys
import time

from helper import parse_dataset_split_chains
from features import get_feature
import logger

logger = logger.get_logger(os.path.basename(__file__))
counter = None

def __compute_feature(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error=False
    errors = []
    try:
        global feature_name
        feat_vals = get_feature(feature_name, input_dir, pdb_id, chain_id)
        output_file = f"{output_dir}/{pdb_id}{chain_id}.txt"
        with open(output_file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in feat_vals))
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as ex:
        error=True
        logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}", exc_info=True)
    finally:
        global counter
        with counter.get_lock():
            idx = counter.value
            counter.value += 1
        if (error):
            errors.append(structure)
            logger.error(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
        else:
            logger.debug(f"{idx}/{total}: {pdb_id} {chain_id} processed")
        return errors

def __pool_init(c):
    ''' store the counter for later use '''
    global counter
    counter = c


dataset_file = ""
output_dir = ""
input_dir = ""
threads = 1
feature_name= ""

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
        feature_name = arg

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
logger.info(f"Computing feature {feature_name} started...")

total = len(dataset)

from multiprocessing import Pool, Value
counter = Value('i', 1)
pool = Pool(int(threads), initializer = __pool_init, initargs = (counter, ))
errors = pool.map(__compute_feature, dataset)
pool.close()
total_errors = [ent for sublist in errors for ent in sublist] #todo delat to lip

if (len(total_errors) == 0):
    logger.info(f"Computing feature {feature_name} finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in total_errors)
    logger.warning(f"Computing feature {feature_name} finished: Some structures were not processed successfully: \n{errors_format}")
logger.debug(f"Finished in {time.time() - start}")