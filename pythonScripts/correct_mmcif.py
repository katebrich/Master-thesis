
import os

dir = '/home/katebrich/Documents/diplomka/PDBe_files/coach420/mmCIF/'

for filename in os.listdir(dir):
    with open(dir + filename, "rt") as fin:
        with open(dir + "temp.txt", "wt") as fout:
            for line in fin:
                fout.write(line.replace('; ', ';\n'))
    os.remove(dir + filename)
    os.rename(dir + "temp.txt", dir + filename)