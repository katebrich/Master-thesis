import getopt
import itertools
import os
import sys
import threading
from helper import eprint, parse_dataset
from features import *


def __compute_feature(str, name_of_feature):
    pdb_id = str[0]
    chain_id = str[1]

    feat_vals = get_feature(name_of_feature, input_dir, pdb_id, chain_id)

    output_file = f"{output_dir}/{pdb_id}{chain_id}.txt"
    with open(output_file, 'w') as f:
        f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in feat_vals))

    with threadLock:
        global counter
        idx = counter
        counter += 1
    print(f"{idx}/{total}: {pdb_id} {chain_id} processed.")


dataset_file = ""
output_dir = ""
input_dir = ""
threads = 1
feature=""

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:i:t:')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
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
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (output_dir == ""):
    eprint("ERROR: Output directory must be specified.")
    sys.exit(1)
if (input_dir == ""):
    eprint("ERROR: Input directory must be specified.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#todo else smazat??

dataset = parse_dataset(dataset_file) #todo co kdyz neni spravny format

total = len(dataset)
counter = 1
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur

print(f"DEBUG: running on {threads} threads.")

if (threads == 1):
    for structure in dataset:
        __compute_feature(feature, structure)
else:
    args = zip(dataset, itertools.repeat(feature, len(dataset)))
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.starmap(__compute_feature, args)
    pool.close()