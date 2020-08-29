import os
from helper import get_fasta_path
from DynaMine.predictor import DynaMine


pdb_id = "2g46"
chain_id = "A"
dataset_name = "test"

data_dir = f"/home/katebrich/Documents/diplomka/TEST_todelete/"
result_dir = f"{data_dir}dynamine/"
fasta_path = f"{data_dir}{pdb_id}{chain_id}.fasta"

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

dynamine = DynaMine(result_dir)
if dynamine.predict(fasta_path, pdb_id, chain_id):
    print(f"DynaMine: {pdb_id} {chain_id} successfully processed.")
else:
    print(f"Error: DynaMine: {pdb_id} {chain_id}")