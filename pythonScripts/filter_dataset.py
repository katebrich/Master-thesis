from parse_dataset import parse_dataset

import os as os
import shutil
from functions import restAPI_get, get_entity_id
import uuid









dataset_name = "joined(mlig)"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.txt'
dataset_dir = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}"

out_files_list = []

dataset = parse_dataset(dataset_file)

i = 0
total = len(dataset)

in_dir_fasta = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/FASTA"
out_dir_fasta = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/FASTA_filtered"  # todo
if os.path.exists(out_dir_fasta):
    shutil.rmtree(out_dir_fasta)
os.makedirs(out_dir_fasta)

in_dir_pdb = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/PDB"
out_dir_pdb = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/PDB_filtered"  # todo
if os.path.exists(out_dir_pdb):
    shutil.rmtree(out_dir_pdb)
os.makedirs(out_dir_pdb)

for structure in dataset:
    i += 1
    pdb_id = structure[0]
    chain_ids = structure[1]

    # filter FASTA file - leave only sequence for this chain
    get_FASTA(in_dir_fasta, out_dir_fasta, pdb_id, chain_ids)


    # filter PDB - leave only these chains + filter relevant ligands
    get_PDB(in_dir_pdb, out_dir_pdb, pdb_id, chain_ids, True)

    print(f"{i}/{total}: {pdb_id} {chain_ids} processed")

dataset_file_dir = os.path.dirname(dataset_file)
with open(f"{dataset_file_dir}/{dataset_name}_prank.ds", 'w') as f:
    for item in out_files_list:
        relpath = os.path.relpath(item, dataset_file_dir)
        f.write(f"{relpath}\n")



