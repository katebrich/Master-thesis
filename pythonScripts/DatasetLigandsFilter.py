import time

import numpy
from Bio.PDB import PDBParser
from scipy.spatial import distance

import Logger
import os

from MOAD import MOAD
from helper import parse_dataset_split_chains, get_pdb_path2, res_mappings_author_to_pdbe, isPartOfChain, \
    getFullAuthorResNum
from multiprocessing import Pool, Value

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class DatasetLigandsFilter:
    filter = ""
    pdb_dir = ""

    def __init__(self, filter):
        self.filter = filter #todo check

    def run(self, orig_dataset_path, filtered_dataset_path, pdb_dir, threads=1):
        self.pdb_dir = pdb_dir # todo check if exists
        start = time.time()
        logger.info(f"Filtering ligands in dataset {orig_dataset_path} with filter '{filter}' started...")

        dataset = parse_dataset_split_chains(orig_dataset_path)

        self.total = len(dataset)

        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        results = pool.map(self.filter_ligands, dataset)
        pool.close()
        total_results = [ent for sublist in results for ent in sublist]

        with open(filtered_dataset_path, 'w') as f:
            f.write('\n'.join('{}\t{}\t{}'.format(x[0], x[1], ','.join(x[2])) for x in results))


        logger.debug(f"Finished in {time.time() - start}")

    def filter_ligands(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        pdb_file = get_pdb_path2(self.pdb_dir, pdb_id, chain_id)
        mappings = dict(
            res_mappings_author_to_pdbe(pdb_id, chain_id))

        ligands_all = []
        AAs = []

        parser = PDBParser(PERMISSIVE=0, QUIET=1)
        structure = parser.get_structure(pdb_id + chain_id, pdb_file)
        chain = structure[0][chain_id]
        # get ligands
        for residue in chain.get_residues():
            if not isPartOfChain(residue, mappings):
                # trim spaces at the beginning of ligand resname; bug in PDB parser
                residue.resname = residue.resname.lstrip()
                ligands_all.append(residue)
            else:
                AAs.append(residue)

        global counter
        with counter.get_lock():
            idx = counter.value
            counter.value += 1
        print(f"{idx}/{self.total}: {pdb_id} {chain_id} processed")

        return (pdb_id, chain_id, self.get_valid_ligands(ligands_all, AAs, pdb_id, chain_id))



    def get_valid_ligands(self, ligands, AAs, pdb_id, chain_id):
        result = []
        skipped = []
        small = []
        center_far = []
        far = []
        if (self.filter == ""):
            for ligand in ligands:
                result.append(ligand.resname)
        elif (self.filter == "p2rank"):
            for ligand in ligands:
                #name of the PDB group is not on the list of ignored groups:
                ignored = ["HOH", "DOD", "WAT", "NAG", "MAN", "UNK", "GLC", "ABA", "MPD", "GOL", "SO4", "PO4"]
                if (ligand.resname in ignored):
                    skipped.append(ligand.resname) #todo debug
                    continue
                #number of ligand atoms is greater than or equal to 5:
                if (len(ligand.child_list) < 5):
                    small.append(ligand.resname)  # todo debug
                    continue
                #distance form the center of the mass of the ligand to the closest protein atom is not greater than 5.5 A #todo tohle pravidlo moc nevychazi
                center = self.__getCenterOfMass(ligand.child_list)
                threshold = 5.5
                success = False
                i = 0
                for AA in AAs:
                    for atom in AA.child_list:
                       if distance.euclidean(atom.coord, center) <= threshold:
                            success = True
                            i += 1
                            break
                if (success == False):
                    center_far.append(ligand.resname)  # todo debug
                    continue
                #distance from any atom of the ligand to the closest protein atom is at least 4 A
                success = False
                for AA in AAs:
                    if (self.__isInDistance(4, AA, ligand)):
                        success = True
                        break
                if (success == False): #no atom in distance 4A was found
                    far.append(ligand.resname)  # todo debug
                    continue

                result.append(ligand.resname)
            print(
                f"{pdb_id} {chain_id}: Remaining {len(result)}/{len(ligands)} - {result}...Center far: {center_far}, Small: {small}, Ignored: {skipped}, Far: {far}") #todo debug only
        elif (self.filter == "MOAD"):
            logger.info(f"Downloading MOAD data file (~18 MB). This may take a while..")
            self.moad = MOAD()
            logger.info(f"MOAD data file downloaded.")


            relevant = self.moad.get_relevant_ligands(pdb_id, chain_id)
            if (relevant == None):
                #raise ValueError("Structure excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5") #todo hezci hlaska, vsechny duvody
              #  print(f"Structure {pdb_id} {chain_id} excluded from MOAD - no valid ligands or not an x-ray structure or has resolution higher than 2.5\n")
                return result # todo
            for ligand in ligands:
                if ((str(ligand.resname), str(ligand.id[1])) in relevant):
                    result.append(ligand.resname)
            #print( f"RELEVANT {pdb_id} {chain_id} : {relevant}")
            if (len(result) != len(relevant)):
                pass
                print(f"ERROR: {pdb_id} {chain_id}:\nRESULT: {result} \nRELEVANT: {relevant}\nLIGANDS: {ligands}\n")
            else:
                pass
               # print(f"OK: {pdb_id} {chain_id}\nRESULT: {result} \nRELEVANT: {relevant}\nLIGANDS: {ligands}\n")

        else:
            raise ValueError(f"Unknown filter level {self.filter}")

        return result

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args

    def __getCenterOfMass(self, atoms):
        totalMass = 0.0
        x = 0
        y = 0
        z = 0
        for a in atoms:
            m = a.mass
            totalMass += m
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
        return (x / totalMass, y / totalMass, z / totalMass)

    def __isInDistance(self, threshold, AA_res, ligand_res):
        for atom1 in AA_res.child_list:
            if atom1.element == 'H': #consider only heavy atoms
                continue
            for atom2 in ligand_res.child_list:
                if atom2.element == 'H':
                    continue
                dist = numpy.linalg.norm(atom1.coord - atom2.coord)
                if dist < threshold:
                    return True
        return False