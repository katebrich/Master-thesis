from parse_dataset import parse_dataset
from Bio import SeqIO
import os as os
from functions import get_fasta_path

#todo parametry
dataset_name = "coach420"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/coach420.ds'
dataset_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420"

dataset = parse_dataset(dataset_file)

i = 1
total = len(dataset)

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    #filter FASTA file - leave only sequence for this chain
    sequence = ""
    file_path = get_fasta_path(dataset_dir, pdb_id, chain_id)
    records = list(SeqIO.parse(file_path, "fasta"))
    for record in records:
        ids = record.description.split('|')
        chains = ids[2].split(' ')
        if chain_id in chains:
            sequence = record.seq
            break

    if sequence == "":
        print(f"Error: sequence not found for {pdb_id} {chain_id}") #todo
        break

    #create temp file
    temp_file_path = file_path + ".temp"
    with open(temp_file_path, 'w') as file:
        file.write(str(sequence.upper()))
    #remove original file
    os.remove(file_path)
    #rename temp file
    os.rename(temp_file_path, file_path)

    print(f"{i}/{total}: {pdb_id} {chain_id} processed")
    i += 1

