import urllib.request
import os, sys
from parse_dataset import parse_dataset
import itertools
from functions import get_pdb_path, get_fasta_path, get_mmcif_path

def download_file(url, output_path):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    with open(output_path, 'wb') as file:
        file.write(responseBody)

def download_structure(str):
    pdb_id = str[0][0]
    total = str[1]
    idx = str[2]
    #chain_id = structure[1]
    # get .pdb file
    url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'
    output_path = os.path.join(output_dir, "PDB", pdb_id + ".pdb") #, chain_id)
    try:
        download_file(url, output_path)
    except: #todo konkretni vyjimky
        print(f"Error: downloading PDB file - {pdb_id}:", sys.exc_info())

    # get .fasta file
    url=f"http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta"
    output_path = os.path.join(output_dir, "FASTA", pdb_id + ".fasta") #chain_id)
    try:
        download_file(url, output_path)
    except:
        print(f"Error: downloading FASTA file - {pdb_id}:", sys.exc_info()[0])

    # get .mmcif file
    #url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif"
    #url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
    #output_path = get_mmcif_path(output_dir, pdb_id, chain_id)
    #download_file(url, output_path)

    #todo vypisovani, ale aby to bylo nejak rychlejsi
   # print(f"{idx}/{total}: {pdb_id} downloaded")


#todo parametry
dataset_name = "coach420"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'

output_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}"
#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)
#if not os.path.exists(f"{output_dir}/FASTA"):
#    os.makedirs(f"{output_dir}/FASTA")
if not os.path.exists(f"{output_dir}/PDB"):
    os.makedirs(f"{output_dir}/PDB")
#if not os.path.exists(f"{output_dir}/mmCIF"):
#    os.makedirs(f"{output_dir}/mmCIF")
    #print(f"Directory {output_dir} created.")

structures = parse_dataset(dataset_file)
idx = range(1, len(structures)+1)
dataset = zip(structures, itertools.repeat(len(structures)), idx)

from multiprocessing.dummy import Pool as ThreadPool
pool = ThreadPool(4)
pool.map(download_structure, dataset)
pool.close()
