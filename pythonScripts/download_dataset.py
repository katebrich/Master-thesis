import urllib.request
import os, sys
from parse_dataset import parse_dataset
from functions import get_pdb_path, get_fasta_path, get_mmcif_path

def download_file(url, output_path):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    with open(output_path, 'wb') as file:
        file.write(responseBody)

#todo parametry
dataset_name = "holo4k_2"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'

output_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(f"{output_dir}/FASTA"):
    os.makedirs(f"{output_dir}/FASTA")
if not os.path.exists(f"{output_dir}/PDB"):
    os.makedirs(f"{output_dir}/PDB")
#if not os.path.exists(f"{output_dir}/mmCIF"):
#    os.makedirs(f"{output_dir}/mmCIF")
    #print(f"Directory {output_dir} created.")

dataset = parse_dataset(dataset_file)

i = 1
total = len(dataset)

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    # get .pdb file
    url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'
    output_path = get_pdb_path(output_dir, pdb_id, chain_id)
    try:
        download_file(url, output_path)
    except: #todo konkretni vyjimky
        print(f"Error: downloading PDB file - {pdb_id} {chain_id}:", sys.exc_info())

    # get .fasta file
    url=f"http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta"
    output_path = get_fasta_path(output_dir, pdb_id, chain_id)
    try:
        download_file(url, output_path)
    except:
        print(f"Error: downloading FASTA file - pdb_id chain_id:", sys.exc_info()[0])

    # get .mmcif file
    #url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif"
    #url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
    #output_path = get_mmcif_path(output_dir, pdb_id, chain_id)
    #download_file(url, output_path)


    # todo upravit soubory, aby tam byl jen ten jeden chain?

    print(f"{i}/{total}: {pdb_id} {chain_id} downloaded")

    i += 1
