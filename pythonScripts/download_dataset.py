import urllib.request
import os
from parse_dataset import parse_dataset
from functions import get_pdb_path, get_fasta_path, get_mmcif_path

def download_file(url, output_path):
    # todo kdyz chyba, vypsat, preskocit!
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()
    with open(output_path, 'wb') as file:
        file.write(responseBody)

#todo parametry
dataset_name = "coach420"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/coach420.ds'

output_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(f"{output_dir}/FASTA"):
    os.makedirs(f"{output_dir}/FASTA")
if not os.path.exists(f"{output_dir}/PDB"):
    os.makedirs(f"{output_dir}/PDB")
if not os.path.exists(f"{output_dir}/mmCIF"):
    os.makedirs(f"{output_dir}/mmCIF")
    #print(f"Directory {output_dir} created.")

dataset = parse_dataset(dataset_file)

i = 1
total = len(dataset)

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    # get .pdb file
    #url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'
    #output_path = get_pdb_path(output_dir, pdb_id, chain_id)
    #download_file(url, output_path)

    # get .fasta file
    #url=f"http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta"
    #output_path = get_fasta_path(output_dir, pdb_id, chain_id)
    #download_file(url, output_path)

    # get .mmcif file
    #url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif"
    url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
    output_path = get_mmcif_path(output_dir, pdb_id, chain_id)
    download_file(url, output_path)

    # todo upravit soubory, aby tam byl jen ten jeden chain?

    print(f"{i}/{total}: {pdb_id} {chain_id} downloaded")

    i += 1
