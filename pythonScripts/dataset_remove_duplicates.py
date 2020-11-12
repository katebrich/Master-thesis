import os
from collections import Counter

from helper import parse_dataset_not_split_chains, parse_dataset_split_chains
import re
from Bio.PDB import *

dataset_name = "holo4k"
datasets_dir = f"/home/katebrich/Documents/diplomka/GitHub/datasets/"
dataset_path = f"{datasets_dir}/{dataset_name}.txt"


dataset = parse_dataset_split_chains(dataset_path)

print(Counter(dataset))

