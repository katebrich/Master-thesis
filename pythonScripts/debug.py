import os

from Bio import SeqIO

from helper import *
from propka.run import single
#from compute_ligand_binding_sites import __compute_ligand_binding_sites

data_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_test"
pdb_id = "1qhi"
chain_id = "A"

def get_FASTA(temp_file, out_dir, pdb_id, chain_ids): #todo zapsat tam i ten header?
    if (chain_ids == '*'): # get all chains #todo tohle tam asi nebude
        url = f"http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta"
        response = restAPI_get(url)
        with open(temp_file, 'wb') as file:
            file.write(response)
        records = list(SeqIO.parse(temp_file, "fasta"))
        for record in records:
            ids = record.description.split('|')
            chains = ids[2].split(' ')
            for chain_id in chains:
                sequence = record.seq
                out_fasta_file_path = os.path.join(out_dir, pdb_id + chain_id + ".fasta")
                with open(out_fasta_file_path, 'w') as file:
                    # todo vypsat i header?
                    file.write(str(sequence.upper()))
    else:
        chains = chain_ids.split(',')
        for chain_id in chains:
            entity_id = get_entity_id(pdb_id, chain_id) #todo mit udelany cache
            url = f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta?entity={entity_id}"
            response = restAPI_get(url)
            data = response.decode().split('\n')
            header = data[0].split('|')
            header[2] = chain_id
            out_fasta_file_path = os.path.join(out_dir, pdb_id + chain_id + ".fasta")
            with open(out_fasta_file_path, 'w') as file:
                file.write(f"{header[0]}|{header[1]}|{header[2]}\n")
                file.write(f"{data[1]}\n")

#__compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_i

get_FASTA(data_dir + "/temp.fasta", data_dir + "/debug", "1g1s", "B")

#print(get_HSE(data_dir, pdb_id, chain_id))