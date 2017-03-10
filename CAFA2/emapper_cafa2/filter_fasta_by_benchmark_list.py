import sys
from Bio import SeqIO
targets = set([line.strip() for line in open('../Supplementary_data/data/benchmark/lists/all.txt')])

for r in SeqIO.parse(open("all_cafa2_seqs.fa"), format="fasta"):
    if r.id in targets:
        print '>%s\n%s' %(r.id, r.seq)
