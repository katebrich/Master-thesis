import gzip
import os

from numpy.core.defchararray import isnumeric

input_dir = "/home/katebrich/Documents/diplomka/P2Rank/datasets_old/joined/conservation/e5i1/scores"
output_dir = "/home/katebrich/Documents/diplomka/GitHub/conservation_results/joined_unchanged"
unchanged_file = "/home/katebrich/Documents/diplomka/GitHub/conservation_results/joined_conservation_unchanged"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

with open(unchanged_file, 'r') as unchanged:
    for unch in unchanged:
        success = False
        for input in os.listdir(input_dir):
            if input.startswith(unch[:5]) or input.startswith(unch[:4] + '_' + unch[4:5]):
                with gzip.open(os.path.join(input_dir, input), 'r') as f:
                    values = []
                    for line in f:
                        line = line.decode('utf-8')
                        columns = line.split()
                        if (len(columns) == 0):
                            continue
                        if (isnumeric(columns[0])):
                            if columns[2][0] != '-':
                                val = float(columns[1])
                                if val == -1000:
                                    val = 0
                                values.append(val)
                    with open(os.path.join(output_dir, f"{unch[:5]}.json"), 'w') as out:
                        out.write('{' + f"\"header\": \">pdb|{unch[:4]}|{unch[4:5]}\", \"database\": \"swissprot\", \"conservation\": {values}" + '}')
                    success = True
                    break
        if not success:
            raise ValueError(f"nothing found for {line}")

