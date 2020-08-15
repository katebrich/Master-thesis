import getopt
import sys
from os.path import isfile
import os

#dataset_file = f'/home/katebrich/Documents/diplomka/datasets/{dataset_name}.txt'
from helper import eprint

dataset_dir = ""
out_path = ""

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:')
except getopt.GetoptError as err:
    eprint(f"ERROR: {err}") #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--data_dir"):
        dataset_dir = arg
    elif opt in ("-o", "--out_path"):
        out_path = arg

if (dataset_dir == ""):
    eprint("ERROR: Dataset directory must be specified.") #todo psat z jakeho skriptu je chyba
    sys.exit(1)
if (out_path == ""):
    eprint("ERROR: Output file pathmust be specified.")
    sys.exit(1)

out_dir = os.path.dirname(dataset_dir) #go one level up
#output_path = os.path.join(out_dir, name)

with open(out_path, 'w') as out_file:
    for f in os.listdir(dataset_dir):
        f_path = os.path.join(dataset_dir, f)
        if isfile(f_path):
            relpath = os.path.relpath(f_path, out_dir)
            out_file.write(f"{relpath}\n")



