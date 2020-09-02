from helper import *
from propka.run import single
from compute_ligand_binding_sites import __compute_ligand_binding_sites

data_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_test"
pdb_id = "1qhi"
chain_id = "A"

__compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_id))

#print(get_HSE(data_dir, pdb_id, chain_id))