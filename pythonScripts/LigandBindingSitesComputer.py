import uuid

from Bio.PDB import *

import SASA
from helper import *
import getopt
import sys
import os
import threading
import logger
import time
from MOAD import MOAD
import numpy

logger = logger.get_logger(os.path.basename(__file__))
counter = None

class LigandBindingSitesComputer():
    dataset_file = ""
    output_dir = ""
    mappings_dir=""
    pdb_dir= ""
    distance_threshold = ""
    SASA_threshold = ""
    filter_level = "p2rank"  # todo parametr
    moad=""
    total = ""

    def __init__(self, dataset_file, output_dir, mappings_dir, PDB_dir, distance_threshold=4, SASA_threshold=0.5):
        self.output_dir = output_dir
        self.mappings_dir = mappings_dir
        self.pdb_dir = PDB_dir
        self.distance_threshold = distance_threshold
        self.SASA_threshold = SASA_threshold
        self.dataset_file = dataset_file

    def run(self, threads):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        dataset = parse_dataset_split_chains(self.dataset_file)

        if (self.filter_level == "MOAD"):
            start = time.time()
            logger.info(f"Downloading MOAD data file (~18 MB)...")
            self.moad = MOAD()
            logger.info(f"MOAD data file downloaded.")
            logger.debug(f"Finished in {time.time() - start}")

        start = time.time()
        logger.info(f"Computing ligand binding sites started...")

        self.total = len(dataset)

        from multiprocessing import Pool, Value
        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.get_ligand_binding_site, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]  # todo delat to lip

        if (len(total_errors) == 0):
            logger.info(f"Computing ligand binding sites finished: All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Computing ligand binding sites finished: Some structures were not processed successfully: \n{errors_format}")
        logger.debug(f"Finished in {time.time() - start}")

    def get_ligand_binding_site(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        errors = []
        error=False
        pdb_path = get_pdb_path2(self.pdb_dir, pdb_id, chain_id)
        try:
            lbs = self.compute_ligand_binding_sites(pdb_id, chain_id, pdb_path)
            output_file = get_lbs_path2(self.output_dir, pdb_id, chain_id)
            with open(output_file, 'w') as f:
                f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in lbs))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error=True
            logger.exception(f"Error while processing {pdb_id} {chain_id}: {ex}", exc_info=True)
        finally:
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append(structure)
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED !")
            else:
                #pass #todo
                logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")
            return errors

    def filter_ligands(self, ligands, pdb_id, chain_id):
        result = []
        skipped = []
        small = []
        center_far = []
        if (self.filter_level == "none"):
            result = ligands
        elif (self.filter_level == "p2rank"):
            for ligand in ligands:
                #name of the PDB group is not on the list of ignored groups:
                ignored = ["HOH", "DOD", "WAT", "NAG", "MAN", "UNK", "GLC", "ABA", "MPD", "GOL", "SO4", "PO4"]
                if (ligand.resname in ignored):
                    skipped.append(ligand.resname) #todo debug
                    continue
                #number of ligand atoms is greater or equal than 5:
                if (len(ligand.child_list) < 5):
                    small.append(ligand.resname)  # todo debug
                    continue
                #distance form the center of the mass of the ligand to the closest protein atom is not greater than 5.5A #todo tohle pravidlo moc nevychazi
                #center = getCenterOfMass(ligand.child_list)
                #threshold = 5.5
                #success = False
                #i = 0
                #for AA in AAs:
                #    #if (success == True):
                #    #    break
                #    for atom in AA.child_list:
                #       if distance.euclidean(atom.coord, center) <= threshold:
                #            success = True
                #            i += 1
                #            break
                #if (success == False):
                #    center_far.append(ligand.resname)  # todo debug
                #    continue
                result.append(ligand)
            #print(
            #    f"{pdb_id} {chain_id}: Remaining {len(result)}/{len(ligands)} - {[x.resname for x in result]}\nFar: {center_far}, Small: {small}, Ignored: {skipped}\n") #todo debug only
        elif (self.filter_level == "MOAD"):
            relevant = self.moad.get_relevant_ligands(pdb_id, chain_id)
            if (relevant == None):
                #raise ValueError("Structure excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5") #todo hezci hlaska, vsechny duvody
              #  print(f"Structure {pdb_id} {chain_id} excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5\n")
                return result # todo
            for ligand in ligands:
                if ((str(ligand.resname), str(ligand.id[1])) in relevant):
                    result.append(ligand)
            #print( f"RELEVANT {pdb_id} {chain_id} : {relevant}")
            if (len(result) != len(relevant)):
                pass
                print(f"ERROR: {pdb_id} {chain_id}:\nRESULT: {result} \nRELEVANT: {relevant}\nLIGANDS: {ligands}\n")
            else:
                pass
               # print(f"OK: {pdb_id} {chain_id}\nRESULT: {result} \nRELEVANT: {relevant}\nLIGANDS: {ligands}\n")

        else:
            raise ValueError(f"Unknown filter level {self.filter_level}")

        return result

    # returns list of tuples:
    # [0] residue number corresponding to PDBe molecule
    # [1] feature value - 0=not binding residue, 1=binding residue
    def compute_ligand_binding_sites(self, pdb_id, chain_id, pdb_file_path):
        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, get_mappings_path2(self.mappings_dir, pdb_id, chain_id)))

        AAs = []
        ligands_all = []

        parser = PDBParser(PERMISSIVE=0, QUIET=1)
        structure = parser.get_structure(pdb_id + chain_id, pdb_file_path)
        chain = structure[0][chain_id]
        #get ligands
        for residue in chain.get_residues():
            if not isPartOfChain(residue, mappings):
                #trim spaces at the beginning of ligand resname; bug in PDB parser
                residue.resname = residue.resname.lstrip() #todo tohle je potreba na MOAD, prehodit jinam, tady smazat?
                ligands_all.append(residue)

        ligands = self.filter_ligands(ligands_all, pdb_id, chain_id) #todo tohle tu asi nebude
        #ligands = ligands_all

        # create temporary .pdb file without HETATM lines to compute Solvent Accessible Surface atoms
        temp_file_path = os.path.join(self.output_dir, f"temp_{uuid.uuid1()}")
        try:
            with open(pdb_file_path, 'r') as orig:
                with open(temp_file_path, 'w') as temp:
                    for line in orig:
                        if line.startswith("ATOM"):
                            temp.write(line)
            parser = PDBParser(PERMISSIVE=0, QUIET=1)
            structure = parser.get_structure(pdb_id + chain_id, temp_file_path)
            chain = structure[0][chain_id]
            sr = SASA.ShrakeRupley()
            sr.compute(structure, level="R")
            for residue in chain.get_residues():
                threshold = 0.05 * sasa_residues[residue.resname]
                if (isPartOfChain(residue, mappings) and residue.sasa > threshold): #leave out residues that are not on the surface
                    AAs.append(residue)
        finally:
            if os.path.exists(temp_file_path):
                os.remove(temp_file_path)

        output = []
        for AA in AAs:
            success = False
            auth_res_num = getFullAuthorResNum(AA.id)
            pdbe_res_num = mappings[auth_res_num]
            for ligand in ligands:
                if (self.isInDistance(AA, ligand)):
                    output.append((pdbe_res_num, 1))
                    success = True
                    break
            if (success == False):
                output.append((pdbe_res_num, 0))
        return output

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args

    #True if there exists a heavy (non-hydrogen) atom on SAS of AA residue which is at maximum distance "threshold" from at least one ligand residue heavy atom
    def isInDistance(self, AA_res, ligand_res):
        # todo optimize, KD tree?
        for atom1 in AA_res.child_list:
            if atom1.element == 'H' or atom1.sasa == 0:  # consider only heavy atoms on Solvent Accessible Surface (SAS)
                continue
            for atom2 in ligand_res.child_list:
                if atom2.element == 'H':
                    continue
                dist = numpy.linalg.norm(atom1.coord - atom2.coord)
                if dist < self.distance_threshold:  # the '-' operator measures distance
                    return True
        return False


sasa_residues = {
    'ALA': 113,   #alanine
    'CYS': 140,   #cysteine
    'ASP': 151,  #aspartic acid
    'GLU': 183,  #glumatic acid
    'PHE': 218,   #phenylalanine
    'GLY': 85,  #glycine
    'HIS': 194,  #histidine
    'ILE': 182,   #isoleucine
    'LYS': 211,  #lysine
    'LEU': 180,   #leucine
    'MET': 204,   #methionine
    'ASN': 158,  #asparagine
    'PRO': 143,  #proline
    'GLN': 189,  #glutamine
    'ARG': 241,  #arginine
    'SER': 122,  #serine
    'THR': 146,  #threonine
    'VAL': 160,   #valine
    'TRP': 259,  #tryptophan
    'TYR': 229,   #tyrosine
}