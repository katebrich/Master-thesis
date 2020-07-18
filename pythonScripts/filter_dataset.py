from parse_dataset import parse_dataset
from Bio import SeqIO
import os as os
from functions import get_fasta_path

#todo parametry
dataset_name = "holo4k"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'
dataset_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}"

dataset = parse_dataset(dataset_file)

i = 0
total = len(dataset)

for structure in dataset:
    i += 1
    pdb_id = structure[0]
    chain_id = structure[1]

    #filter FASTA file - leave only sequence for this chain
    sequence = ""
    file_path = get_fasta_path(dataset_dir, pdb_id, chain_id)
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
    except Exception as e:
        print(f"Error: {pdb_id} {chain_id}")
        print(str(e))
        continue

    for record in records:
        ids = record.description.split('|')
        chains = ids[2].split(' ')
        if chain_id in chains:
            sequence = record.seq
            break

    if sequence == "":
        print(f"Error: sequence not found for {pdb_id} {chain_id}") #todo
        continue

    #create temp file
    temp_file_path = file_path + ".temp"
    with open(temp_file_path, 'w') as file:
        file.write(str(sequence.upper()))
    #remove original file
    os.remove(file_path)
    #rename temp file
    os.rename(temp_file_path, file_path)

    print(f"{i}/{total}: {pdb_id} {chain_id} processed")


