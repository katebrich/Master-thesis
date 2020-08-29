from Bio.PDB import *
from helper import *
import getopt
import sys
import os
import threading
import logger
import time

logger = logger.get_logger(os.path.basename(__file__))

def __get_ligand_binding_site(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error=False
    pdb_path = get_pdb_path(data_dir, pdb_id, chain_id)
    try:
        lbs = __compute_ligand_binding_sites(pdb_id, chain_id, pdb_path)
        output_file = get_lbs_path(data_dir, pdb_id, chain_id)
        with open(output_file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in lbs))
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

# returns list of tuples:
# [0] residue number corresponding to PDBe molecule
# [1] feature value - 0=not binding residue, 1=binding residue
def __compute_ligand_binding_sites(pdb_id, chain_id, pdb_file_path):
    parser = PDBParser(PERMISSIVE=0, QUIET=1) #todo
    # parser = MMCIFParser()
    structure = parser.get_structure(pdb_id + chain_id, pdb_file_path)
    chain = structure[0][chain_id]

    mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))

    AAs = []
    ligands = []

    for residue in chain.get_residues():
        if (residue.id[0] == ' '):  # hetero flag is empty
            AAs.append(residue)
        elif (residue.id[0] != 'W'):
            # todo filter relevant ligands
            # ligand_atoms = ligand_atoms + residue.child_list
            ligands.append(residue)

    output = []
    for AA in AAs:
        success = False
        auth_res_num = str(AA.id[1])
        if (AA.id[2] != ' '):
            auth_res_num += str(AA.id[2])  # insertion code
        pdbe_res_num = mappings[auth_res_num]
        for ligand in ligands:
            if (isInDistance(THRESHOLD, AA, ligand)):
                output.append((pdbe_res_num, 1))
                success = True
                break;
        if (success == False):
            output.append((pdbe_res_num, 0))

    return output


dataset_file = ""
data_dir = ""
threads = 1
THRESHOLD = 4  # in angstroms #todo parametr

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

lbs_dir = os.path.join(data_dir, "lbs") #todo helper?

if not os.path.exists(lbs_dir):
    os.makedirs(lbs_dir)

dataset = parse_dataset_split_chains(dataset_file) #todo co kdyz neni spravny format

start = time.time()
logger.info(f"Computing ligand binding sites started...")

total = len(dataset)
counter = 1
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur
errors = []

if (threads == 1):
    for structure in dataset:
        __get_ligand_binding_site(structure)
else:
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.map(__get_ligand_binding_site, dataset)
    pool.close()

if (len(errors) == 0):
    logger.info(f"Computing ligand binding sites finished: All structures processed successfully.")
else:
    errors_format = '\n'.join('%s %s' % x for x in errors)
    logger.warning(f"Computing ligand binding sites finished: Some structures were not processed successfully: \n{errors_format}")
logger.debug(f"Finished in {time.time() - start}")