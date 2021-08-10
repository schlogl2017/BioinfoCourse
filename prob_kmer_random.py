from random import choice
import numpy as np

def randseq(n):
    return ''.join([choice(['A', 'T', 'C', 'G']) for i in range(n)])

counts = []
i = 0
while len(counts) < 5 or (3 * (np.std(counts) / np.sqrt(i))) > 0.0001:
    i += 1
    kmer = randseq(9)
    kcount = 0
    for j in range(500):
        seq = randseq(1000)
        if kmer in seq:
            kcount += seq.count(kmer)
    counts.append(kcount)
    print(i, kcount, np.mean(counts), (3 * (np.std(counts) / np.sqrt(i))))
