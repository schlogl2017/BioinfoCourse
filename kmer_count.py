#!usr/bin/env python
import os
import argparse




k = 6
fastaFile = '/path/to/some/seqs.fa'
kmerCmd = 'kmer-counter --fasta --k=%d %s' % (k, fastaFile)
try:
    output = subprocess.check_output(kmerCmd, shell=True)
    result = {}
    for line in output.splitlines():
        (header, counts) = line.strip().split('\t')
        header = header[1:]
        kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in counts.split(' ')])
        result[header] = kmers
    sys.stdout.write("%s" % (str(result)))
except subprocess.CalledProcessError as error:
    sys.stderr.write("%s" % (str(error)))
