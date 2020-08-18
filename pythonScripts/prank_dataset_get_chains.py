import os
from helper import parse_prank_dataset
import re
from Bio.PDB import *

dataset_name = "mine_holo4k_part2"
datasets_dir = f"/home/katebrich/Documents/diplomka/P2Rank_with_csv_feature/datasets_old/"
dataset_path = f"{datasets_dir}/{dataset_name}.ds"

output_dir = "/home/katebrich/Documents/diplomka/datasets/from_prank/"
output_path = f"{output_dir}/{dataset_name}.txt"
#create the output directory or remove its contents, if it already exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

dataset = parse_prank_dataset(dataset_path)
i = 1
total = len(dataset)

with open(output_path, 'w') as file:
    for line in dataset:
        #print(f"Processing {line}")
        pdb_ids = re.findall(r"[0-9][A-Za-z0-9]{3}", str(line))
        if (len(pdb_ids) != 1):
            print(f"Error: unable to determine PDB ID from the file name: {line}")
            continue
        pdb_id = pdb_ids[0]
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        structure = parser.get_structure(pdb_id, os.path.join(datasets_dir, line))
        ids = []
        for chain in structure.get_chains():
            chain_id = chain.id
            if (chain_id == " "): #remove chains without name
                continue
            non_ligand = False

            # filter out chians that content only ligands:
            for res in chain.get_residues():
                if (res.id[0] == ' '):  # hetero flag is empty
                    non_ligand = True
            if non_ligand:
                ids.append(chain_id)

        if len(ids) == 0:
            ids.append('A') #default chain

        file.write(f"{pdb_id}\t{ ','.join(ids)}\n")

        print(f"{i}/{total}: structure {line} processed.")
        i += 1