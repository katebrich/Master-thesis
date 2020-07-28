from parse_dataset import parse_dataset
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import os as os
import shutil
from functions import get_fasta_path

#todo parametry
dataset_name = "joined(mlig)"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.txt'
#dataset_dir = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}"

out_files_list = []

class ChainSelect(Select):
    def __init__(self, chain):
        super(Select, self).__init__()
        self.chain = chain
    def accept_residue(self, residue):
        if residue.parent.id ==self.chain and residue.id[0] != 'W':
            return 1
        else:
            return 0
    def accept_model(self, model):
        if model.id == 0:
            return 1
        else:
            return 0

def filter_PDB(in_dir, out_dir, pdb_id, chain_ids, filter_ligands=True):
    in_pdb_file_path = os.path.join(in_dir, pdb_id + ".pdb") #todo windows?
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, in_pdb_file_path)
    chains_list = []
    if (chain_ids == '*'):
        for chain in structure.get_chains():
            chains_list.append(chain.id)
    else:
        chains_list = chain_ids.split(',')
    #todo model
    for chain_id in chains_list: #todo chain_id = *
        out_pdb_file_path = os.path.join(out_dir, pdb_id + chain_id + ".pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_pdb_file_path, ChainSelect(chain_id), preserve_atom_numbering=True)
        out_files_list.append(out_pdb_file_path)

def filter_FASTA(in_dir, out_dir, pdb_id, chain_ids):
    in_fasta_file_path = os.path.join(in_dir, pdb_id + ".fasta") #todo windows?
    try:
        records = list(SeqIO.parse(in_fasta_file_path, "fasta"))
    except Exception as e:
        print(f"Error: {pdb_id} {chain_ids}")
        print(str(e))
        return

    #todo
    #if (chain_ids == '*'):
    #   for chain in structure.get_chains():
    #        chains_list.append(chain.id)
    #else:
    chains_list = chain_ids.split(',')

    for record in records:
        ids = record.description.split('|')
        chains = ids[2].split(' ')
        for chain_id in chains:
            if chain_id in chains_list:
                sequence = record.seq
                chains_list.remove(chain_id)
                out_fasta_file_path = os.path.join(out_dir, pdb_id + chain_id + ".fasta")
                with open(out_fasta_file_path, 'w') as file:
                    #todo vypsat i header?
                    file.write(str(sequence.upper()))
    if len(chains_list) != 0:
        print(f"Error: sequence not found for {pdb_id} {chains_list}")  # todo

    ## create temp file
    #temp_file_path = file_path + ".temp"
    #with open(temp_file_path, 'w') as file:
    #    file.write(str(sequence.upper()))
    ## remove original file
    #os.remove(file_path)
    ## rename temp file
    #os.rename(temp_file_path, file_path)


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
    filter_FASTA(in_dir_fasta, out_dir_fasta, pdb_id, chain_ids)


    # filter PDB - leave only these chains + filter relevant ligands
    filter_PDB(in_dir_pdb, out_dir_pdb, pdb_id, chain_ids, True)

    print(f"{i}/{total}: {pdb_id} {chain_ids} processed")

dataset_file_dir = os.path.dirname(dataset_file)
with open(f"{dataset_file_dir}/{dataset_name}_prank.ds", 'w') as f:
    for item in out_files_list:
        relpath = os.path.relpath(item, dataset_file_dir)
        f.write(f"{relpath}\n")



