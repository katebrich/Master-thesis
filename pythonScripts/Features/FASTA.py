from Bio import SeqIO
from helper import *
from Features.AA_properties import *

def get_AA_properties(scores_dict, data_dir, pdb_id, chain_id):
    fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
    seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
    return get_AA_scores(scores_dict, seq)

class AA():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq) + 1):
            AA = seq[i - 1]
            if AA in ['B','J','O','U','X','Z']:
                continue
            feature_vals.append((i, AA))
        return feature_vals

class AA_pairs():
    def get_values(self, data_dir, pdb_id, chain_id):
        fasta_file = get_fasta_path(data_dir, pdb_id, chain_id)
        seq = list(SeqIO.parse(fasta_file, "fasta"))[0]
        feature_vals = []
        for i in range(1, len(seq)):
            aa_1 = seq[i - 1]
            aa_2 = seq[i]
            feature_vals.append((i, aa_1 + aa_2))
        return feature_vals

class Hydropathy():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(hydropathy_kyte_doolitle, data_dir, pdb_id, chain_id)

class Polarity():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(polarity, data_dir, pdb_id, chain_id)

class PolarityBinary():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(polarity_binary, data_dir, pdb_id, chain_id)

class Charged():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(charged, data_dir, pdb_id, chain_id)

class Aromaticity():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(aromaticity, data_dir, pdb_id, chain_id)

class MolecularWeight():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(molecular_weight, data_dir, pdb_id, chain_id)

class HBondAtoms():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(H_bond_atoms, data_dir, pdb_id, chain_id)

class AARatio():
    def get_values(self, data_dir, pdb_id, chain_id):
        return get_AA_properties(aa_ratio, data_dir, pdb_id, chain_id)