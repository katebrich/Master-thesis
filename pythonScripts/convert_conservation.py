import gzip
import os

from numpy.core.defchararray import isnumeric

input_dir = "/home/katebrich/Documents/diplomka/P2Rank/datasets_old/coach420/conservation/e5i1/scores"
output_dir = "/home/katebrich/Documents/diplomka/P2Rank/datasets_old/coach420/conservation/e5i1/values_only"

for input in os.listdir(input_dir):
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
        print(values)
