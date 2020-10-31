import os

from Bio import SeqIO

import helper

from helper import *
from propka.run import single
#from compute_ligand_binding_sites import __compute_ligand_binding_sites

#get_entity_id("1a4k", "H")

#__compute_ligand_binding_sites(pdb_id, chain_id, get_pdb_path(data_dir, pdb_id, chain_i

#get_FASTA(data_dir + "/temp.fasta", data_dir + "/debug", "1g1s", "B")

#print(get_HSE(data_dir, pdb_id, chain_id))

class_name="test.Testik"

def my_import(name):
    components = name.split('.')
    mod = __import__(components[0])
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

klass = my_import(class_name)
some_object = klass()

some_object.baf("Kacka")