import sys
import getopt
from Bio.PDB import PDBParser
from helper import parse_prank_dataset
import os
import pandas as pd
import re
from features import get_feature
from helper import res_mappings_author_to_pdbe
import csv
import shutil

#try:
#    opts, args = getopt.getopt(sys.argv[1:], 'd:f:')
#except getopt.GetoptError as err:
#    # print help information and exit:
#    print(err) # will print something like "option -a not recognized"
#    #todo print help
#    sys.exit(2)
#for opt, arg in opts:
#    if opt in ("-h", "--help"):
#        # todo print help
#       sys.exit()
#    elif opt in ("-d", "--directory"):
#        dir = arg
#    elif opt in ("-f", "--file"):
#        file = arg
#    #todo else unknown option

label = "prank_hydropathy_without_fix"
dataset_name = "chen11"
dataset_path = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}_prank.ds'
data_dir = f"/home/katebrich/Documents/diplomka/datasets/{dataset_name}/" #todo parametr
dataset_file_dir = os.path.dirname(dataset_path)
output_dir = os.path.join(data_dir + label)
#create the output directory or remove its contents, if it already exists
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir)

dataset = parse_prank_dataset(dataset_path)
features = ["hydropathy"] #todo
defaults = [0] #todo

i = 1
total = len(dataset)
#aa_codes = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"] #TODO smazat az se odstrani chyba v pranku

for line in dataset:
    print(f"Processing {line}")
    df = pd.DataFrame(columns=["chain", "ins. code", "seq. code"] + features)
    pdb_ids = re.findall(r"[0-9][A-Za-z0-9]{3}", str(line)) #todo test
    if (len(pdb_ids) != 1):
        print(f"Error: unable to determine PDB ID from the file name: {line}")
        continue
    pdb_id = pdb_ids[0]
    parser = PDBParser(PERMISSIVE=0, QUIET=1)  # todo
    structure = parser.get_structure(pdb_id, os.path.join(dataset_file_dir, line))
    for chain in structure.get_chains():
        chain_id = chain.id
        feature_vals = []
        for feature in features:
            feature_vals.append(dict(get_feature(feature, data_dir, pdb_id, chain_id)))
        mappings = res_mappings_author_to_pdbe(pdb_id, chain_id)
        for residue in chain.get_residues():
            if False:#(residue.id[0][2:] in aa_codes): #todo smazat fix pro prank
                if (residue.id[2].isspace()):
                    ins_code = ""
                else:
                    ins_code = residue.id[2]
                seq_code = residue.id[1]
                res_num = str(seq_code)+str(ins_code)
                feat_tuple = tuple(defaults)
                df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code, seq_code) + feat_tuple
            elif (residue.id[0].isspace() or residue.id[0] == "H_MSE" ): #skip hetero-residues #todo selenomethionine???
                if (residue.id[2].isspace()): #biopython returns a space instead of empty string
                    ins_code = ""
                else:
                    ins_code = residue.id[2]
                seq_code = residue.id[1]
                res_num = mappings[str(seq_code)+str(ins_code)]
                feat_list = []
                #todo missing feature vals for res_num
                for j in range(0, len(features)):
                    val = feature_vals[j].get(res_num, defaults[j])
                    feat_list.append(val)
                feat_tuple = tuple(feat_list)
               # print(chain_id, ins_code, seq_code, feat_tuple)
                df.loc[0 if pd.isnull(df.index.max()) else df.index.max() + 1] = (chain_id, ins_code, seq_code) + feat_tuple
    output_path = os.path.join(output_dir, os.path.basename(line) + ".csv")
    with open(output_path, 'w') as file:
        file.write(df.to_csv(index=False, quoting=csv.QUOTE_ALL))

    print(f"{i}/{total}: structure {line} processed.")
    i += 1


