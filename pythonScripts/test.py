import os

from Bio import SeqIO

dir1 = "/home/katebrich/Documents/diplomka/P2Rank/datasets/coach420_10_27/FASTA"
dir2 = "/home/katebrich/Documents/diplomka/P2Rank/datasets_old/coach420/fasta"

unchanged = []
changed = []

for filename1 in os.listdir(dir1):
    for filename2 in os.listdir(dir2):
        if filename2.startswith(filename1[:5]):
            seq1 = list(SeqIO.parse(os.path.join(dir1, filename1), "fasta"))[0]
            seq2 = list(SeqIO.parse(os.path.join(dir2, filename2), "fasta"))[0]
            if seq1.seq == seq2.seq:
                unchanged.append(filename1[:5])
            else:
                changed.append(filename1[:5])
            break


print('\n'.join('%s.fasta' % x for x in changed))



