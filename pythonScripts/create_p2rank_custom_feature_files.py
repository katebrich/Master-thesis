import sys
import getopt
import threading
import traceback
from Bio.PDB import PDBParser
from helper import parse_prank_dataset, eprint, parse_dataset_split_chains, get_pdb_path
import os
import pandas as pd
import re
from features import default_values
from helper import res_mappings_author_to_pdbe
import csv
import shutil
import numpy as np

def read_feature_vals(feature, data_dir, pdb_id, chain_id):
    path = os.path.join(data_dir, "features", feature, f"{pdb_id}{chain_id}.txt")
    feature_vals = np.genfromtxt(path, delimiter=' ')
    return dict(feature_vals)

def create_file(structure):
    pdb_id = structure[0]
    chain_id = structure[1]
    error = False

    try:
        df = pd.DataFrame(columns=["chain", "ins. code", "seq. code"] + features)
        parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
        pdb_path = get_pdb_path(input_dir, pdb_id, chain_id)
        struct = parser.get_structure(pdb_id, pdb_path)

        feature_vals = []
        for feature in features:
            feature_vals.append(read_feature_vals(feature, input_dir, pdb_id, chain_id))  # todo

        mappings = dict(res_mappings_author_to_pdbe(pdb_id, chain_id, os.path.join(mappings_cache_dir, f"{pdb_id}{chain_id}.txt"))) #todo cache
        #print(mappings)
        for residue in struct.get_residues():
            if (residue.id[0][2:] in aa_codes):  # todo smazat fix pro prank
                if (residue.id[2].isspace()):
                    ins_code = ""
                else:
                    ins_code = residue.id[2]
                seq_code = residue.id[1]
                #res_num = str(seq_code) + str(ins_code)
                feat_tuple = tuple(defaults)
                # print(chain_id, ins_code, seq_code, feat_tuple)
                df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code,
                                                                                  seq_code) + feat_tuple
            elif (residue.id[0].isspace() or residue.id[0] == "H_MSE"):  # skip hetero-residues #todo selenomethionine???
                if (residue.id[2].isspace()):  # biopython returns a space instead of empty string
                    ins_code = ""
                else:
                    ins_code = residue.id[2]
                seq_code = residue.id[1]
                res_num = mappings[str(seq_code) + str(ins_code)]
                feat_list = []
                # todo missing feature vals for res_num
                for j in range(0, len(features)):
                    val = feature_vals[j].get(res_num, defaults[j])
                    feat_list.append(val)
                feat_tuple = tuple(feat_list)
                # print(chain_id, ins_code, seq_code, feat_tuple)
                df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code,
                                                                                  seq_code) + feat_tuple
        output_path = os.path.join(output_dir, os.path.basename(pdb_path) + ".csv")
        with open(output_path, 'w') as file:
            file.write(df.to_csv(index=False, quoting=csv.QUOTE_ALL))
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as ex:
        error=True
        print(f"ERROR: downloading {pdb_id} {chain_id}: {ex}")  #todo test jestli se dostane do chyboveho vystupu shell skriptu
        traceback.print_exception(type(ex), ex, ex.__traceback__)
    finally:
        with threadLock:
            global counter
            idx = counter
            counter += 1
        if (error):
            print(f"{idx}/{total}: {pdb_id} {chain_id} NOT PROCESSED !")
        else:
            print(f"{idx}/{total}: {pdb_id} {chain_id} processed")

#todo zrychlit to!!!

dataset_file = ""
output_dir = ""
input_dir = ""
threads = 1 #todo
features={}

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:i:t:')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-i", "--input_dir"): #folder with PDB files
        input_dir = arg
    elif opt in ("-t", "--threads"):
        threads = arg #todo check if threads >= 1, int
    elif opt in ("-f", "--features"):
        features = arg.split(',') #todo check

if (dataset_file == ""):
    eprint("ERROR: Dataset must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (output_dir == ""):
    eprint("ERROR: Output directory must be specified.")
    sys.exit(1)
if (input_dir == ""):
    eprint("ERROR: Input directory must be specified.")
    sys.exit(1)
#todo check features

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#todo else smazat??

dataset = parse_dataset_split_chains(dataset_file) #todo co kdyz neni spravny format

defaults = [default_values[f] for f in features]
print("DEFAULTS:", defaults)

mappings_cache_dir = os.path.join(input_dir, "mappings")
aa_codes = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"] #TODO smazat az se odstrani chyba v pranku
total = len(dataset)
counter = 1
threadLock = threading.Lock() #todo otestovat o kolik to bude rychlejsi bez toho locku a vypisovani processed struktur

print(f"DEBUG: running on {threads} threads.")

if (threads == 1):
    for structure in dataset:
        create_file(structure)
else:
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(int(threads))
    pool.map(create_file, dataset)
    pool.close()



