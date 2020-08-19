import getopt
import os
import sys
import threading
import traceback

from helper import eprint, parse_dataset_split_chains, res_mappings_author_to_pdbe
import logger

logger = logger.get_logger(os.path.basename(__file__))

def create_mappings_cache(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error = False
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
        with threadLock:
            global counter
            idx = counter
            counter += 1
        if (error):
            errors.append(structure)
            logger.error(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
            #todo zapsat nekam chybu
        else:
            logger.debug(f"{idx}/{total}: {pdb_id} {chain_id} processed")

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
    logger.error("Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (output_dir == ""):
    logger.error("Output directory must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


dataset = parse_dataset_split_chains(dataset_file)  #todo co kdyz neni spravny format

logger.info(f"Creating mapping cache started...")

total = len(dataset)
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur
counter = 1
errors = []

if (threads == 1):
    for structure in dataset:
        create_mappings_cache(structure)
else:
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.map(create_mappings_cache, dataset)
    pool.close()

if (len(errors) == 0):
    logger.info(f"Creating mapping cache finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in errors)
    logger.warning(f"Creating mapping cache finished: Some structures were not processed successfully: \n {errors_format}")