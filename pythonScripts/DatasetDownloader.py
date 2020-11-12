import os, sys

from MOAD import MOAD
from helper import parse_dataset_not_split_chains
import threading
import traceback
import getopt
from helper import restAPI_get, get_entity_id
import uuid
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import Logger
import time
from multiprocessing import Pool, Value

logger = Logger.get_logger(os.path.basename(__file__))
counter = None

class DatasetDownloader():
    dataset_file = ""
    output_dir = ""

    total=""
    output_PDB=""
    output_FASTA=""

    #ligands_filter:
    def __init__(self, dataset_file, output_dir):
        self.dataset_file = dataset_file
        self.output_dir = output_dir


    def run(self, threads=1):
        #self.threads = threads
        self.output_PDB = os.path.join(self.output_dir, "PDB")
        self.output_FASTA = os.path.join(self.output_dir, "FASTA")

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.output_PDB):
            os.makedirs(self.output_PDB)
        if not os.path.exists(self.output_FASTA):
            os.makedirs(self.output_FASTA)

        dataset = parse_dataset_not_split_chains(self.dataset_file)  # todo co kdyz neni spravny format

        start = time.time()
        logger.info(f"Downloading structures to {self.output_dir} started...")

        self.total = len(dataset)

        counter = Value('i', 1)
        pool = Pool(int(threads), initializer=self.__pool_init, initargs=(counter,))
        errors = pool.map(self.download_structure, dataset)
        pool.close()
        total_errors = [ent for sublist in errors for ent in sublist]

        if (len(total_errors) == 0):
            logger.info(f"Downloading structures finished: All structures downloaded successfully.")  #todo cas rovnou sem
        else:
            errors_format = '\n'.join('%s %s' % x for x in total_errors)
            logger.warning(
                f"Downloading structures finished: Some structures were not downloaded successfully: \n{errors_format}")
        logger.debug(f"Finished in {time.time() - start}")

    def get_PDB(self, temp_file, out_dir, pdb_id, chain_ids):
        url = f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent'
        response = restAPI_get(url)
        with open(temp_file, 'wb') as file:
            file.write(response)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        structure = parser.get_structure(pdb_id, temp_file)
        chains_list = chain_ids.split(',')
        for chain_id in chains_list:
            out_pdb_file_path = os.path.join(out_dir, pdb_id + chain_id + ".pdb")
            io = PDBIO()
            io.set_structure(structure)
            io.save(out_pdb_file_path, PDBSelector(chain_id), preserve_atom_numbering=True)

    def get_FASTA(self, out_dir, pdb_id, chain_ids):  # todo zapsat tam i ten header?
        chains = chain_ids.split(',')
        for chain_id in chains:
            entity_id = get_entity_id(pdb_id, chain_id)
            url = f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta?entity={entity_id}"
            response = restAPI_get(url)
            data = response.decode().split('\n')
            header = data[0].split('|')
            header[2] = chain_id
            out_fasta_file_path = os.path.join(out_dir, pdb_id + chain_id + ".fasta")
            with open(out_fasta_file_path, 'w') as file:
                file.write(f"{header[0]}|{header[1]}|{header[2]}\n")
                file.write(f"{data[1]}\n")

    def download_structure(self, structure):
        pdb_id = structure[0]
        chain_ids = structure[1]
        error = False
        errors = []
        temp_file = os.path.join(self.output_PDB, f"temp_{uuid.uuid1()}")
        try:
            # temp_file = os.path.join(output_FASTA, f"temp_{uuid.uuid1()}")
            self.get_FASTA(self.output_FASTA, pdb_id, chain_ids)
            # os.remove(temp_file)
            self.get_PDB(temp_file, self.output_PDB, pdb_id, chain_ids)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.debug(f"Error while downloading {pdb_id} {chain_ids}: {ex}", exc_info=True)
        finally:
            if os.path.exists(temp_file):
                os.remove(temp_file)
            global counter
            with counter.get_lock():
                idx = counter.value
                counter.value += 1
            if (error):
                errors.append(structure)
                logger.error(f"{idx}/{self.total}: {pdb_id} {chain_ids} NOT PROCESSED ! See log for more details.")
            else:
                logger.debug(f"{idx}/{self.total}: {pdb_id} {chain_ids} processed")
            return errors

    def __pool_init(self, args):
        ''' store the counter for later use '''
        global counter
        counter = args



class PDBSelector(Select): #todo tohle popsat v praci
    chain = ""
    def __init__(self, chain):
        super(Select, self).__init__()
        self.chain = chain
    def accept_residue(self, residue):
        if residue.parent.id == self.chain and residue.id[0] != 'W': #do not save water hetatoms
            return 1
        else:
            return 0
    def accept_atom(self, atom):
        if atom.element == 'H' or atom.element == 'D': #do not save hydrogens (or deuterium)
            return 0
        else:
            return 1
    def accept_model(self, model):
        if model.id == 0:
            return 1
        else:
            return 0
