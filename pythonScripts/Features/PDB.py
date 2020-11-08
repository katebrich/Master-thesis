from Bio.PDB import PDBParser, HSExposure
from helper import *


class HSE_up():
    def get_values(self, data_dir, pdb_id, chain_id):
        pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]
        HSExposure.HSExposureCB(model)
        feature_vals = []
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
        for r in model.get_residues():
            if not (r.id[0].isspace()):  # is HETATM
                continue
            try:
                hse = int(r.xtra["EXP_HSE_B_U"]) #todo proc vychazi stejne Up i Down?
                auth_res_num = getFullAuthorResNum(r.id)
                pdbe_res_num = mappings[auth_res_num]
                feature_vals.append((pdbe_res_num, hse))
            except:
                pass  # some residues are incomplete in PDB file (i.e. CB atom is missing) and no HSE value is returned

        return feature_vals

class HSE_down():
    def get_values(self, data_dir, pdb_id, chain_id):
        pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]
        HSExposure.HSExposureCB(model)
        feature_vals = []
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
        for r in model.get_residues():
            if not (r.id[0].isspace()):  # is HETATM
                continue
            try:
                hse = int(r.xtra["EXP_HSE_B_D"])
                auth_res_num = getFullAuthorResNum(r.id)
                pdbe_res_num = mappings[auth_res_num]
                feature_vals.append((pdbe_res_num, hse))
            except:
                pass  # some residues are incomplete in PDB file (i.e. CB atom is missing) and no HSE value is returned

        return feature_vals

class ExposureCN():
    def get_values(self, data_dir, pdb_id, chain_id):
        pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]
        HSExposure.ExposureCN(model, radius=10.0)  # radius as in P2Rank protrusion feature
        feature_vals = []
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
        for r in model.get_residues():
            if not (r.id[0].isspace()):  # is HETATM
                continue
            try:
                cn = int(r.xtra["EXP_CN"])
                auth_res_num = getFullAuthorResNum(r.id)
                pdbe_res_num = mappings[auth_res_num]
                feature_vals.append((pdbe_res_num, cn))
            except:
                raise ValueError(
                    f"Error CN: {pdb_id} {chain_id}: residue {r}")  # todo jestli se to pro cely dataset nestane, tak asi muzu smazat try-except

        return feature_vals

class BFactor():
    def get_values(self, data_dir, pdb_id, chain_id):
        pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)
        structure = parser.get_structure(pdb_id, pdb_file)  # todo udelat funkci get_chain?
        chain = structure[0][chain_id]
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
        feature_vals = []
        for residue in chain.get_residues():
            if not isPartOfChain(residue, mappings):
                continue
            sum = 0
            for atom in residue.child_list:
                sum += atom.bfactor
            b_factor = round(sum / len(residue.child_list), 3)
            auth_res_num = getFullAuthorResNum(residue.id)
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, b_factor))
        return feature_vals

class BFactorCA():
    def get_values(self, data_dir, pdb_id, chain_id):
        pdb_file = get_pdb_path(data_dir, pdb_id, chain_id)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)
        structure = parser.get_structure(pdb_id, pdb_file)  # todo udelat funkci get_chain?
        chain = structure[0][chain_id]
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path(data_dir, pdb_id, chain_id)))
        feature_vals = []
        for residue in chain.get_residues():
            if not isPartOfChain(residue, mappings):
                continue
            bfactor = None
            for atom in residue.child_list:
                if atom.id == 'CA':
                    if (bfactor != None):
                        raise ValueError(
                            f"Error: more C alpha atoms in {pdb_id} {chain_id} residue {residue.id}")  # todo only debug
                    bfactor = atom.bfactor
                    # todo zrychlit, kdyz najdu C alpha -> break
            if bfactor == None:
                print(f"Error: no C alpha in {pdb_id} {chain_id} residue {residue.id}")  # todo only debug
                continue
            auth_res_num = getFullAuthorResNum(residue.id)
            pdbe_res_num = mappings[auth_res_num]
            feature_vals.append((pdbe_res_num, bfactor))
        return feature_vals