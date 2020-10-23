import getopt
import os
import sys
import threading
import traceback
import time

from helper import eprint, parse_dataset_split_chains, res_mappings_author_to_pdbe
import logger

logger = logger.get_logger(os.path.basename(__file__))
counter = None

def create_mappings_cache(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error = False
    errors = []
    try:
        mappings = res_mappings_author_to_pdbe(pdb_id, chain_id)
        output_file = f"{output_dir}/{pdb_id}{chain_id}.txt"
        with open(output_file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in mappings))
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
            logger.error(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
            #todo zapsat nekam chybu
        else:
            logger.debug(f"{idx}/{total}: {pdb_id} {chain_id} processed")
        return errors

def __pool_init(args):
    ''' store the counter for later use '''
    global counter
    counter = args

dataset_file = ""
output_dir = ""
threads = 1

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:t:')
except getopt.GetoptError as err:
    logger.error(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-t", "--threads"): #todo check if threads >= 1, int
        threads = arg

if (dataset_file == ""):
    logger.error("Dataset must be specified.")
    sys.exit(1)
if (output_dir == ""):
    logger.error("Output directory must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


dataset = parse_dataset_split_chains(dataset_file)

start = time.time()
logger.info(f"Creating mapping cache started...")

total = len(dataset)

from multiprocessing import Pool, Value
counter = Value('i', 1)
pool = Pool(int(threads), initializer=__pool_init, initargs=(counter,))
errors = pool.map(create_mappings_cache, dataset)
pool.close()
total_errors = [ent for sublist in errors for ent in sublist]

if (len(total_errors) == 0):
    logger.info(f"Creating mapping cache finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in total_errors)
    logger.warning(f"Creating mapping cache finished: Some structures were not processed successfully: \n{errors_format}")
logger.debug(f"Finished in {time.time() - start}")