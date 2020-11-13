import uuid

from Bio.PDB import *
from helper import *
import getopt
import sys
import os
import threading
import Logger
import time
from scipy.spatial import distance
import SASA
from shutil import copyfile

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

def __get_sasa(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    errors = []
    error=False
    pdb_path = get_pdb_path(data_dir, pdb_id, chain_id)
    try:
        # create temporary .pdb file without HETATM lines
        temp_file_path = os.path.join(sasa_dir, f"temp_{uuid.uuid1()}")
        with open(pdb_path, 'r') as orig:
            with open(temp_file_path, 'w') as temp:
                for line in orig:
                    if line.startswith("ATOM"):
                        temp.write(line)

        sasa = __compute_sasa_residues(pdb_id, chain_id, temp_file_path)
        output_file = get_sasa_path(data_dir, pdb_id, chain_id)
        with open(output_file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in sasa))
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
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)
        return errors

# returns list of tuples:
# [0] residue number corresponding to PDBe molecule
# [1] feature value - 0=not solvent accessible residue, 1=solvent accessible residue
def __compute_sasa_residues(pdb_id, chain_id, pdb_file_path):
    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))

    #parse structure without ligands
    parser = PDBParser(PERMISSIVE=0, QUIET=1)
    structure = parser.get_structure(pdb_id + chain_id, pdb_file_path)
    chain = structure[0][chain_id]

    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R")

    output = []
    for residue in chain.get_residues():
        auth_res_num = getFullAuthorResNum(residue.id)
        pdbe_res_num = mappings[auth_res_num]
        output.append((pdbe_res_num, residue.sasa))

    return output


def __pool_init(args):
    ''' store the counter for later use '''
    global counter
    counter = args


dataset_file = ""
data_dir = ""
threads = 1
RADIUS = 1.4


#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:i:t:')
except getopt.GetoptError as err:
    logger.error(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-i", "--data_dir"):
        data_dir = arg
    elif opt in ("-t", "--threads"):
        threads = arg #todo check if threads >= 1, int

if (dataset_file == ""):
    logger.error("Dataset must be specified.")
    sys.exit(1)
if (data_dir == ""):
    logger.error("Input directory must be specified.")
    sys.exit(1)

sasa_dir = os.path.join(data_dir, "sasa") #todo helper?

if not os.path.exists(sasa_dir):
    os.makedirs(sasa_dir)

dataset = parse_dataset(dataset_file) #todo co kdyz neni spravny format


start = time.time()
logger.info(f"Computing solvent accessible surface residues started...")

total = len(dataset)

from multiprocessing import Pool, Value
counter = Value('i', 1)
pool = Pool(int(threads), initializer = __pool_init, initargs = (counter, ))
errors = pool.map(__get_sasa, dataset)
pool.close()
total_errors = [ent for sublist in errors for ent in sublist] #todo delat to lip

if (len(total_errors) == 0):
    logger.info(f"Computing solvent accessible surface residues finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in total_errors)
    logger.warning(f"Computing solvent accessible surface residues finished: Some structures were not processed successfully: \n{errors_format}")
logger.debug(f"Finished in {time.time() - start}")