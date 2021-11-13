# quick script to count the numbers of kmers across multiple genomes. Saves results for all in outFile

import sys
import os
import glob
from kmerAnalysis import doKmerCount

if len(sys.argv) != 5:
	sys.exit("USAGE: kmerCountWindows.py fastaMap k outFile")

# sys.argv[1] = Data/Genomes_splitted
fastaMap = glob.glob(f'{sys.argv[1]}/*/{sys.argv[2]}/*_fna.gz')
k = int(sys.argv[3])
outFile = sys.argv[4]

# read fasta information

with open(outFile, 'w') as of:
	first = True
	for fname in fastaMap:
	    name = fname.split('/')[2]
		data = doKmerCount(fname, k)
		if first:
			# write header of file - names of tetras
			of.write(','.join(data.keys())+'\n')
		# write kmer counts 
		of.write(name+',' + ','.join([str(x) for x in data.values()]) + '\n')
		first = False
