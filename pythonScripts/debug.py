import os
import Config
from Bio import SeqIO

import helper

from helper import *
from propka.run import single
#from compute_ligand_binding_sites import __compute_ligand_binding_sites

#get_entity_id("1a4k", "H")

#__compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_i

#get_FASTA(data_dir + "/temp.fasta", data_dir + "/debug", "1g1s", "B")

#print(get_HSE(data_dir, pdb_id, chain_id))

path = Config.get_feature_path("unp_PTM")
from pydoc import locate
feature_class = locate(path)
#feature_class = class_import(path)
feature_class.get_values()
      #  feat_vals = feature_class.get_values(input_dir, pdb_id, chain_id)