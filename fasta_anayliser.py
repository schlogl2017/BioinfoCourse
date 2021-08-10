#!usr/bin/env python
# -*- coding: utf-8 -*-

class fasta(object):
    def __init__(self, filename):
    self.filename = filename
    self.num_sequences = None
    self.sequences = {} #{seq_id} = sequence

def parse_file(self):
    **"""Reads in the sequence data contained in the filename associated with this  instance of the class.
    Stores both the sequence and the optional comment for each ID."""**
    with open(self.filename) as f:
        return f.read().split('>')[1:]
def get_info(self):
    **"""Returns a description of the class in a pretty string format. 
    The description should include the filename for the instance and the number of sequences."""**
    for line in file(self.filename, 'r'):
        if line.startswith('>'):
            self.num_sequences += 1
    return self.num_sequences
def compute_gc_content(self,some_id):
    **"""compute the gc conent for sequence ID some_id. If some_id, return an appropriate error values"""**
    baseFrequency = {}
    for line in file(self.filename, 'r'):
        if not line.startswith(">"):
        for base in sequence:
            baseFrequency[base] = baseFrequency.get(base,0)+1
            items = baseFrequency.items()
            items.sort()
        for i in items:
            gc=(baseFrequency['G'] + baseFrequency['C'])/float(len(sequence))
    return gc

def sequence_statistics(self):
      **"""returns a dictionary containing
         The average sequence length
         The average gc content"""**
    baseFrequency = {}
    for line in file(self.filename, 'r'):
        if not line.startswith(">"):
        for base in sequence:
            baseFrequency[base] = baseFrequency.get(base,0)+1
            items = baseFrequency.items()
            items.sort()
        for i in items:
            gc=(baseFrequency['G'] + baseFrequency['C'])/float(len(sequence))
            aveseq=sum(len(sequence))/float(self.count)
    return (gc,aveseq)


def get_all_kmers(self, k=8):
    **"""get all kmer counts of size k in fasta file. Returns a dictionary with keys equal to the kmers
    and values equal to the counts"""
    t={}
    for x in range (len(self.sequence)+1-k):
        kmer=self.sequence[x:x+k]
        t[kmer]=f.get(kmer,0)+1
        kmers = get_all_kmers(k=8)
    return(t,len(kmers))

def query_sequence_id(self, some_id):
    **"""query sequence ids for some_id. If some_id does not exist in the class, return
    a string error message"""**
    for line in file(self.filename, 'r'):
    if id in line:
        print "The id exists"
    else:
        print "The id does not exist"