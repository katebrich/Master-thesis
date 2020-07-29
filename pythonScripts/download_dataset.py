import urllib.request
import os, sys
from parse_dataset import parse_dataset
import itertools
from functions import get_pdb_path, get_fasta_path
import threading
import traceback
import shutil
import getopt
from functions import eprint, restAPI_get


def download_file(url, output_path):
    response = restAPI_get(url)
    with open(output_path, 'wb') as file:
        file.write(response)

def download_structure(str):
    try:
        pdb_id = str[0]
        chain_ids = str[1]

        # get .pdb file
        url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'
        output_path = os.path.join(output_dir, "PDB", pdb_id + ".pdb") #, chain_id)
        try:
            download_file(url, output_path)
        except:
            eprint(f"ERROR: downloading PDB file - {pdb_id}") #todo test jestli se dostane do chyboveho vystupu shell skriptu
            raise

        # get .fasta file
        url=f"http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta"
        output_path = os.path.join(output_dir, "FASTA", pdb_id + ".fasta") #chain_id)
        try:
            download_file(url, output_path)
        except: #t
            eprint(f"ERROR: downloading FASTA file - {pdb_id}")

        #todo vypisovani, ale aby to bylo nejak rychlejsi
        with threadLock:
            global global_counter
            idx = global_counter
            global_counter += 1

        print(f"{idx}/{total}: {pdb_id} downloaded")
    except:
        print() #todo test
        # todo log error structures


dataset_file = ""
output_dir = ""
ligands_filter = None

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:l')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-l", "--ligands_filter"):
        ligands_filter = arg #todo check possible values
#todo parameter number of threads

if (dataset_file == ""):
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (output_dir == ""):
    eprint("ERROR: Output directory must be specified.")
    sys.exit(1)

#remove output directory and its contents, if it exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.makedirs(f"{output_dir}/FASTA")
os.makedirs(f"{output_dir}/PDB")

structures = parse_dataset(dataset_file)

# global variables
total = len(structures)     #todo otestovat jestli to multithreading zrychluje
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur
global_counter = 1

from multiprocessing.dummy import Pool as ThreadPool
pool = ThreadPool(4)
pool.map(download_structure, structures)
pool.close()

#todo vypsat kolik bylo chyb a kolik bylo stazeno uspesne