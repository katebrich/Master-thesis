from Bio.PDB import *
import array as array
from functions import isInDistance

def get_ligand_binding_sites(pdb_id, chain_id, pdb_file_path):
    parser = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure(pdb_id + chain_id, pdb_file_path)
    chain = structure[0][chain_id]

    AAs = []
    #ligand_atoms = []
    ligands = []

    for residue in chain.get_residues():
        if (residue.id[0] == ' '): #hetero flag is empty
            AAs.append(residue)
        elif (residue.id[0] != 'W'):
            #todo filter relevant ligands
            #ligand_atoms = ligand_atoms + residue.child_list
            ligands.append(residue)

    THRESHOLD = 4 # in angstroms

    output = array.array('B', [0]*len(AAs))
    count = 0

    for AA in AAs:
        for ligand in ligands:
            if (isInDistance(THRESHOLD, AA, ligand)):
                output[count] = 1
                break;
        count += 1

    return output

