import urllib.request
import os
from parse_dataset import parse_dataset

dataset_name = "coach420"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/coach420.ds'

output_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Directory {output_dir} created.")

dataset = parse_dataset(dataset_file)

i = 1
total = len(dataset)

for structure in dataset:
    pdb_id = structure[0]
    chain_id = structure[1]

    url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'

    #todo kdyz chyba, preskocit!
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        responseBody = f.read()

    with open(f'{output_dir}/{pdb_id}{chain_id}.pdb', 'wb') as file:
        file.write(responseBody)

    print(f"{i}/{total}: {pdb_id} {chain_id} downloaded")

    i += 1
