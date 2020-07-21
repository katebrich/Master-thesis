from parse_dataset import parse_dataset
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import os as os
from functions import get_fasta_path

class ChainSelect(Select):
    def __init__(self, chain):
        super(Select, self).__init__()
        self.chain = chain
    def accept_residue(self, residue):
        if residue.parent.id ==self.chain and residue.id[0] != 'W':
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

    for chain_id in chains_list: #todo chain_id = *
        out_pdb_file_path = os.path.join(out_dir, pdb_id + chain_id + ".pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_pdb_file_path, ChainSelect(chain_id), preserve_atom_numbering=True)


def filter_FASTA(pdb_id, chain_id):
    sequence = ""
    file_path = get_fasta_path(dataset_dir, pdb_id, chain_id)
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
    except Exception as e:
        print(f"Error: {pdb_id} {chain_id}")
        print(str(e))
        return

    for record in records:
        ids = record.description.split('|')
        chains = ids[2].split(' ')
        if chain_id in chains:
            sequence = record.seq
            break

    if sequence == "":
        print(f"Error: sequence not found for {pdb_id} {chain_id}")  # todo
        return

    # create temp file
    temp_file_path = file_path + ".temp"
    with open(temp_file_path, 'w') as file:
        file.write(str(sequence.upper()))
    # remove original file
    os.remove(file_path)
    # rename temp file
    os.rename(temp_file_path, file_path)


#todo parametry
dataset_name = "coach420"
dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.ds'
dataset_dir = f"/home/katebrich/Documents/diplomka/PDBe_files/{dataset_name}"

in_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420/PDB"
out_dir = "/home/katebrich/Documents/diplomka/PDBe_files/coach420/PDB_filtered"  # todo
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dataset = parse_dataset(dataset_file)

i = 0
total = len(dataset)

for structure in dataset:
    i += 1
    pdb_id = structure[0]
    chain_ids = structure[1]

    #filter FASTA file - leave only sequence for this chain
    #filter_FASTA(pdb_id, chain_id)

    #filter PDB - leave only these chains + filter relevant ligands
    filter_PDB(in_dir, out_dir, pdb_id, chain_ids, True)

    print(f"{i}/{total}: {pdb_id} {chain_ids} processed")



