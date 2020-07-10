from Bio.PDB import *
import array as array

from Bio.PDB.MMCIF2Dict import *

from functions import isInDistance

#returns list of tuples:
    # [0] residue number corresponding to PDBe molecule
    # [1] feature value - 0=not binding residue, 1=binding residue
def get_ligand_binding_sites(pdb_id, chain_id, mmcif_file_path):
    #parser = PDBParser(PERMISSIVE=0)
    #pdb_id = "1a4k"
    #chain_id = 'H'
    #mmcif_file_path = f"/home/katebrich/Downloads/1a4k_updated.cif"
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id + chain_id, mmcif_file_path)
    chain = structure[0][chain_id]

    #mappings from author residue number to pdbe molecule residue number
    mmcif_dict = MMCIF2Dict(mmcif_file_path)
    group_PDB = mmcif_dict['_atom_site.group_PDB']
    pdbe_seq_id = mmcif_dict['_atom_site.pdbe_label_seq_id']
    author_seq_id = mmcif_dict['_atom_site.auth_seq_id']
    author_ins_code = mmcif_dict['_atom_site.pdbx_PDB_ins_code']
    author_chain_id = mmcif_dict['_atom_site.auth_asym_id']

    #todo smazat, jen kontrola
    if (len(pdbe_seq_id) != len(group_PDB) or len(pdbe_seq_id) != len(author_seq_id) or len(pdbe_seq_id) != len(author_ins_code) or len(pdbe_seq_id) != len(author_chain_id) ):
        print(f"Error: pdbe_seq_id not same lenght as author_seq_id")
        return

    #todo nekam si ukladat a priste cist ze souboru, nedelat pokazde znovu?
    mappings = []
    for i in range(0, len(pdbe_seq_id)):
        if (group_PDB[i] != "ATOM" or author_chain_id[i] != chain_id):
            continue;
        val = pdbe_seq_id[i]
        key = author_seq_id[i]
        if (author_ins_code[i] != '?'):
            key += author_ins_code[i]
        mappings.append((key, val))
    mappings_dict = dict(mappings)

    AAs = []
    ligands = []

    for residue in chain.get_residues():
        if (residue.id[0] == ' '): #hetero flag is empty
            AAs.append(residue)
        elif (residue.id[0] != 'W'):
            #todo filter relevant ligands
            #ligand_atoms = ligand_atoms + residue.child_list
            ligands.append(residue)

    THRESHOLD = 3.5 # in angstroms

    output = []
    for AA in AAs:
        success = False
        key = str(AA.id[1])
        if (AA.id[2] != ' '):
            key += str(AA.id[2]) #insertion code
        pdbe_res_num = mappings_dict[key]
        for ligand in ligands:
            if (isInDistance(THRESHOLD, AA, ligand)):

                output.append((pdbe_res_num, 1))
                success = True
                break;
        if (success == False):
            output.append((pdbe_res_num, 0))

    return output

