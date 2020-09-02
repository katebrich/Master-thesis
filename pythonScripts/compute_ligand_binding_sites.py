from Bio.PDB import *
from helper import *
import getopt
import sys
import os
import threading
import logger
import time
from MOAD import MOAD
from scipy.spatial import distance

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

def filter_ligands(ligands, AAs, pdb_id, chain_id):
    result = []
    skipped = []
    small = []
    center_far = []
    if (filter_level == "none"):
        result = ligands
    elif (filter_level == "p2rank"):
        for ligand in ligands:
            #name of the PDB group is not on the list of ignored groups:
            ignored = ["HOH", "DOD", "WAT", "NAG", "MAN", "UNK", "GLC", "ABA", "MPD", "GOL", "SO4", "PO4"]
            if (ligand.resname in ignored):
                skipped.append(ligand.resname) #todo debug
                continue
            #number of ligand atoms is greater or equal than 5:
            if (len(ligand.child_list) < 5):
                small.append(ligand.resname)  # todo debug
                continue
            #distance form the center of the mass of the ligand to the closest protein atom is not greater than 5.5A #todo tohle pravidlo moc nevychazi
            center = getCenterOfMass(ligand.child_list)
            threshold = 5.5
            success = False
            i = 0
            for AA in AAs:
                #if (success == True):
                #    break
                for atom in AA.child_list:
                    if distance.euclidean(atom.coord, center) <= threshold:
                        success = True
                        i += 1
                        break
            if (success == False):
                center_far.append(ligand.resname)  # todo debug
                continue
            result.append(ligand)
            print(
                f"{pdb_id} {chain_id}: Remaining {len(result)}/{len(ligands)} - {[x.resname for x in result]}\nFar: {center_far}, Small: {small}, Ignored: {skipped}\n") #todo debug only
    elif (filter_level == "MOAD"):
        relevant = MOAD.get_relevant_ligands(pdb_id, chain_id)
        if (relevant == None):
            raise ValueError("Structure excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5") #todo hezci hlaska, vsechny duvody
        for ligand in ligands:
            if ((str(ligand.resname), str(ligand.id[1])) in relevant):
                result.append(ligand)
        print( f"RELEVANT {pdb_id} {chain_id} : {relevant}")

    else:
        raise ValueError(f"Unknown filter level {filter_level}")

    return result

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
    ligands_all = []
    ids = [] #todo only for debug

    for residue in chain.get_residues():
        if (residue.id[0] == ' '):  # hetero flag is empty
            AAs.append(residue)
        elif (is_aa(residue, standard=False)): #HETATM, but nonstandard AA code (MSE, LYZ etc.)
            #check if mapping exist for this residue number:
            auth_res_num = getFullAuthorResNum(residue.id)
            #auth_res_num = str(residue.id[1])
            #if (residue.id[2] != ' '):
            #    auth_res_num += str(residue.id[2])  # insertion code
            if (auth_res_num in mappings):
                AAs.append(residue) #is part of the chain, not ligand
            else:
                ligands_all.append(residue) #is truly a ligand
        else:
            # todo filter relevant ligands
            ligands_all.append(residue)
            #ids.append(residue.id[0])

    ligands = filter_ligands(ligands_all, AAs, pdb_id, chain_id)

    #print(f"LIGANDS: {pdb_id} {chain_id}: {ids}")
    output = []
    for AA in AAs:
        success = False
        auth_res_num = getFullAuthorResNum(AA.id)
        #auth_res_num = str(AA.id[1])
        #if (AA.id[2] != ' '):
        #    auth_res_num += str(AA.id[2])  # insertion code
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
filter_level = "p2rank"

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

#todo debug
#dataset_file="/home/katebrich/Documents/diplomka/datasets/chen11.txt"
#data_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_chen11"
#threads = 1
#filter_level="MOAD"

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

if (filter_level == "MOAD"):
    start = time.time()
    logger.info(f"Downloading MOAD data file (~18 MB)...")
    MOAD = MOAD()
    logger.info(f"MOAD data file downloaded.")
    logger.debug(f"Finished in {time.time() - start}")

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