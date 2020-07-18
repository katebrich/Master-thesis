import os
from functions import get_fasta_path
from predictor import DynaMine


pdb_id = "1a2k"
chain_id = "C"
dataset_name = "coach420"

data_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}/"
result_dir = f"{data_dir}dynamine/"
fasta_path = f"{data_dir}FASTA_original/{pdb_id}{chain_id}.fasta"

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

dynamine = DynaMine(result_dir)
if dynamine.predict(fasta_path):
    print(f"DynaMine: {pdb_id} {chain_id} successfully processed.")
else:
    print(f"Error: DynaMine: {pdb_id} {chain_id}")