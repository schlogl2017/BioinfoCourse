#!/usr/bin/env python3

import sys, math
from itertools import product
import fastaReader as fred

########################################################################
# File: analyzeKmerFreq.py
#  executable: analyzeKmerFreq.py
# Purpose: Find underrepresented 3-9mers in a FASTA file via Z-score
#   stderr: errors and status
#   stdout: list of 3-9mers found in genome sorted by Z-score within kmer size
#          
# Author: Christopher Lee
#               
########################################################################

class AnalyzeKmers():
    '''
    Compute Kmer Z-scores.

    Attributes:
        kmerDict - Dictionary of all observed kmer sequences in input genome.

    Output:
        A tab-sep list of kmers, their observed counts, expected count and Z-score
        in that order.
    '''
    def __init__(self, minSize, maxSize):
        self.observedKmers = {}
        self.expectedKmers = {}
        self.minSize = minSize
        self.maxSize = maxSize
        self.validNucs = 'ATCG'
        self.genomeSize = 0
    
    def findKmers(self,dnaSeq):
        '''Iterates through the sequence and finds all kmers between min and max size.'''
        self.genomeSize += len(dnaSeq)
        for i in range(0, self.genomeSize): # test all bases so we don't miss overlapping kmers
            for j in range(self.minSize-2, self.maxSize+1): # find kmers of size min-2 : max for later
                if i+j <= self.genomeSize:
                    for char in dnaSeq[i:i+j]: # don't count non canonical bases
                        if char not in self.validNucs:
                            continue
                    if dnaSeq[i:i+j] in self.observedKmers.keys(): 
                        self.observedKmers[dnaSeq[i:i+j]] += 1
                    else:
                        self.observedKmers[dnaSeq[i:i+j]] = 2 # add 2 to match example output, should really be a 1 though
   
    def permuteKmers(self):
        '''Finds list of all kmers of size min to max, mostly a wrapper around itertools product lib'''
        for k in range(self.minSize,self.maxSize+1):
            kmers = list(''.join(x) for x in product(self.validNucs, repeat=k))
            self.expectedKmers.update( {kmer:0 for kmer in kmers} )

    def computeExpectedCounts(self, kmer):
        '''Finds the expected count for a specifc kmer via checking the observed counts of the 
           sub-kmers.'''
        left = kmer[:-1]
        right = kmer[1:]
        middle = kmer[1:-1]
        try:
            expCount = (self.observedKmers[left] * self.observedKmers[right]) / (
                        self.observedKmers[middle])#/ self.genomeSize)
        except KeyError:
            expCount = 0
        return expCount

    def computeZScore(self, cutoff=0):
        '''Returns a list of lists:
               [[kmer, observed count, expected count, z-score above cutoff]]
        
           Args: 
               cutoff: Defaults to zero since we are looking for under-represented seqs
                       with negative z-scores
        '''
        outList = []
        kmersToAnalyze = [kmer for kmer in self.expectedKmers.keys()] # if (
        for kmer in kmersToAnalyze:
            expCount = self.computeExpectedCounts(kmer)
            prob = expCount / self.genomeSize
            try:
                obsCount = self.observedKmers[kmer]
            except KeyError:
                obsCount = 0 
            mean  = self.genomeSize*prob
            stdDev = math.sqrt((self.genomeSize*prob*(1-prob)))
            try:
                zScore = ((obsCount - mean) / stdDev)
            except ZeroDivisionError:
                zScore = 0 
            if zScore < cutoff:
                outList.append([kmer, obsCount, expCount, zScore])
        return outList

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    attributes:
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Analyze Kmer Frequency by Z-score', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options]  < input > output' 
                                             )
        self.parser.add_argument('-c', '--cutoff', type=float, action='store', help='The Z-score cutoff, must be a float')
        self.parser.add_argument('-m', '--maxMotif', type=int, default=8, action='store', help='The maximum Kmer size to evaluate')
        self.parser.add_argument('-l', '--minMotif', type=int, default=3, action='store', help='The minimum Kmer size to evaluate')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
  

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg 


def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    kmerAnalyzer = AnalyzeKmers(myCommandLine.args.minMotif, myCommandLine.args.maxMotif)
    faReader = fred.FastaReader()
    for header,seq in faReader.readFasta():
        kmerAnalyzer.findKmers(seq)
    kmerAnalyzer.permuteKmers()
    out = kmerAnalyzer.computeZScore()
    sortedOut = sorted(sorted(out,key=lambda item: item[3]), key=lambda item: len(item[0]), reverse=True)
    for kmer in range(len(sortedOut)):
        print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(sortedOut[kmer][0], sortedOut[kmer][1], sortedOut[kmer][2], sortedOut[kmer][3]))

if __name__ == "__main__":
    main()
