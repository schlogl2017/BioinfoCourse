#!usr/bin/env python
# script to know (for each aligned position, the percentage of 
# amioacid species, eg, 80% A, 10% G etc


import sys
from Bio import AlignIO
from collections import Counter


alignment = AlignIO.read(sys.argv[1], 'phylip')

for i in range(alignment.get_alignment_length()):
    print(Counter(aln[:, i))
