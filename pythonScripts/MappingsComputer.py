import getopt
import os
import sys
import threading
import traceback
import time

from helper import parse_dataset_split_chains, res_mappings_author_to_pdbe, getStructuresFromDirectory
import Logger

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class MappingsComputer:
    output_dir = ""
    dataset_file = ""
    total=""

    def __init__(self, dataset_file, output_dir):
        self.output_dir = output_dir
        self.dataset_file = dataset_file

    def run(self, threads):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset_split_chains(self.dataset_file)

        start = time.time()
        logger.info(f"Creating mapping cache started...")

        self.total = len(dataset)

        from multiprocessing import Pool, Value
        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.create_mappings_cache, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Creating mapping cache finished: All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Creating mapping cache finished: Some structures were not processed successfully: \n{errors_format}")
        logger.debug(f"Finished in {time.time() - start}")

    def create_mappings_cache(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        errors = []
        try:
            mappings = res_mappings_author_to_pdbe(pdb_id, chain_id)
            output_file = f"{self.output_dir}/{pdb_id}{chain_id}.txt"
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
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED !")
                #todo zapsat nekam chybu
            else:
                logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args


'''
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
'''

