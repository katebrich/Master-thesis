import getopt
import os
import sys
import threading
import traceback

from helper import eprint, parse_dataset_split_chains, res_mappings_author_to_pdbe


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
        print(f"ERROR: processing {pdb_id} {chain_id}: {ex}")
        traceback.print_exception(type(ex), ex, ex.__traceback__)
    finally:
        with threadLock:
            global counter
            idx = counter
            counter += 1
        if (error):
            print(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
            #todo zapsat nekam chybu
        else:
            print(f"{idx}/{total}: {pdb_id} {chain_id} processed")

dataset_file = ""
output_dir = ""
threads = 1

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:t:')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-t", "--threads"): #todo check if threads >= 1, int
        threads = arg

if (dataset_file == ""):
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (output_dir == ""):
    eprint("ERROR: Output directory must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


dataset = parse_dataset_split_chains(dataset_file)  #todo co kdyz neni spravny format

total = len(dataset)
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur
counter = 1

if (threads == 1):
    for structure in dataset:
        create_mappings_cache(structure)
else:
    print(f"DEBUG: running on {threads} threads.")
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.map(create_mappings_cache, dataset)
    pool.close()