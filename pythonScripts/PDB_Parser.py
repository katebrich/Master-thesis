from Bio.PDB import *

path = '/home/katebrich/Documents/diplomka/PDBe_files/pdb1igy.ent'
chain_id = "B"

parser = PDBParser(PERMISSIVE=0)
structure = parser.get_structure('str', path)
chain = structure[0][chain_id]
residues = list(chain.get_residues())

print(residues[0])
print(residues[123])
print(residues[125])

for residue in residues:
    print(residue)
    print(residue.get_resname())  # return the residue name (eg. 'GLY')
    print(residue.is_disordered())  # 1 if the residue has disordered atoms
    #print(residue.get_segid())  # return the SEGID
    #residue.has_id(name)\

    #measure distances of atoms
    #at0 = residue['CA']
    #at1 = residue['N']
    #distance = at0 - at1
    #print(distance)

    #test if AA
    print(is_aa(residue))

