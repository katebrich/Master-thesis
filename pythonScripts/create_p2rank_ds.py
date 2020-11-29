import getopt
import sys
from os.path import isfile
import os

from Bio.PDB import PDBParser

import Logger
from helper import isPartOfChain

'''
def get_ligands_codes(path):
    # get all ligands from PDB file
    ligands=[]
    parser = PDBParser(PERMISSIVE=0, QUIET=1)
    structure = parser.get_structure("structure", path)
    for chain in structure[0]:
        for residue in chain.get_residues():
            if not residue.id[0].isspace(): #todo or residue.id[0] == "H_MSE"): #hetatm
                # trim spaces at the beginning of ligand resname; bug in PDB parser
                residue.resname = residue.resname.lstrip()
                ligands.append(residue)
    return ligands
'''

logger = Logger.get_logger(os.path.basename(__file__))

pdb_dir = ""
out_path = ""
list_ligands=False
dataset_path=""

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'p:o:d:l')
except getopt.GetoptError as err:
    logger.error(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-p", "--pdb_dir"):
        pdb_dir = arg
    elif opt in ("-o", "--out_path"):
        out_path = arg
    elif opt in ("-d", "--dataset_path"):
        dataset_path = arg
    elif opt in ("-l", "--list_ligands"): #todo bez parametru
        list_ligands = True

#todo
if (pdb_dir == ""):
    logger.error("Dataset directory must be specified.")
    sys.exit(1)
if (out_path == ""):
    logger.error("Output file path must be specified.")
    sys.exit(1)

out_dir = os.path.dirname(out_path) #go one level up
#output_path = os.path.join(out_dir, name)

ligands_dict={}
if list_ligands:
    with open(dataset_path, 'r') as f:
        for line in f:
            columns = line.split()
            ligands_dict[columns[0]+columns[1]] = columns[2]

with open(out_path, 'w') as out_file:
    if list_ligands:
        out_file.write(f"HEADER: protein ligand_codes\n\n")
    for f in os.listdir(pdb_dir):
        f_path = os.path.join(pdb_dir, f)
        if isfile(f_path):
            relpath = os.path.relpath(f_path, out_dir)
            if list_ligands:
                ligands = ligands_dict[f[:5]]
                out_file.write(f"{relpath}\t{ligands}\n")
            else:
                out_file.write(f"{relpath}\n")



