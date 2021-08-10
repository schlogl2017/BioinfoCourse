#!usr/bin/env python


import gzip
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from collections import defaultdict
from itertools import groupby


def is_header(line):
    """
    Function to check if a line is a header in a fasta file.
    
    Inputs:
    
        line - a line from a fasta file.
    
    Outputs:
    
        Returns True if the line is a header in a fasta file, other wise False.
    """
    return line[0] == '>'


def parse_fasta(filename):
    """
    Function to read a fasta file.
    
    Inputs:
    
        filename - a path to a fasta file.
    
    Outputs:
    
        Yields the header and the sequence from a fasta file.
    """
    # checks if the file is compressed or not
    if filename.endswith('.gz'):
        # when compressed
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')
    # read the fasta file
    with opener(filename) as f:
        # group the header and the sequence
        fasta_iter = (it[1] for it in groupby(f, is_header))
        # gets the headers and the sequence
        # then yields both
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequences = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequences.upper()


def get_name_sequence_biopython(filename):
    """
    Function to read a fasta file using the biopython module.
    
    Inputs:
    
        filename - a path to a fasta file.
    
    Outputs:
    
        Yields the header and the sequence from a fasta file.
    """
    # reads the file
    for rec in SeqIO.parse(gzip.open(filename, "rt"), "fasta"):
        # gets the name/id 
        name = rec.id
        # gets the sequence and yields both
        seq = str(rec.seq)
        yield name, seq


def fasta_item_counter(filename):
    """
    Function to read and count the sequences in a fasta file.
    
    Inputs:
    
        filename - a path to a fasta file.
    
    Outputs:
    
        count - integer representing the number of sequences in a fasta file.
    """
    # checks if the file is compressed or not
    # counts the number of sequences and returns the count
    if filename.endswith('.gz'):
        count = sum(g for g, _ in groupby(gzip.open(filename, 'rt'), key=is_header))
    else:
        count = sum(g for g, _ in groupby(open(filename, 'r'), key=is_header))
    return count


def str_punctuation_strip(word):
    """
    Function to clean up a word from any kind of punctuation.
    
    Inputs:
    
        word - a string representing a word.
    
    Outputs:
    
        word - a string representing a cleaned word.
    """
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    # iterates through the word characters
    for _ in word:
        # iterates through the punctuation characters
        for p in punctuation:
            # replace the puntiation for a spaces
            word = word.replace(p, ' ')
    # returns the cleaned word
    return ''.join(word.strip().split())


def get_fasta_headers(filename):
    """
    Function to read a fasta file and returning the headers.
    
    Inputs:
    
        filename - a path to a fasta file.
    
    Outputs:
    
        headers - if is a multi fasta file (more than one sequence).
        header - if a one sequence fasta file.
    """
    # initialize the array to receive the headers
    headers = []
    # read the headers
    for name, _ in parse_fasta(filename):
        # add the header to the array
        headers.append(name)
    # if multi fasta
    if len(headers) > 1:
        return headers
    # if single fasta
    else:
        return name


def split_sequences_from_fasta_file(filename):
    """
    Function that get a multifasta file and split it in to single fasta files sequence.
    Input:
    
        filename - a full path tho the multi fasta file. Ex; a assembly or a reference
                   ncbi file of a complete bacterial genome, with chromosome and
                   plasmid sequences.
    
    Output:
    
        plasmids, chromosome - a tuple of dictionary-like objects with chromosome and
                               plasmids splitted from a input multi fasta file.                   
    """
    plasmids = defaultdict(list)
    chromosome = defaultdict(list)
    # get sequence ID from the filename. Ex. 
    # ['Data/Assemblies/Acidiphilium/GCF_000202835.1.fna.gz']
    # GCF_000202835.1
    id_seq = filename.split('/')[-1][:15]
    # read the file and get the headers and sequences
    for header, seq in parse_fasta(filename):
        # checking for plasmid sequences in the file
        # and get them in the map
        if 'plasmid' in header:
            ids = id_seq + '_plas'
            plasmids[ids] = plasmids.get(ids, [])
            plasmids[ids].append(seq)
        else:
            # checking for the chromossome sequence and mapping
            ids = id_seq + '_chr'
            chromosome[ids] = chromosome.get(ids, [])
            chromosome[ids].append(seq)
    # join the plasmid sequences
    for name, seq in plasmids.items():
        plasmids[name] = ''.join(seq)
    # join the chromossome sequences
    for name, seq in chromosome.items():
        chromosome[name] = ''.join(seq)
    return plasmids, chromosome


def write_fasta_file(sequence_dict, out_file, wrap=80):
    """Write sequences to a fasta file.

    Inputs:
    
        sequence_dict - dictionary-like object with the sequences indexed by sequence id.
        out_file - path and name of the fasta file to write the sequences to.
        wrap - a number of characters by line (default = 80).
    
    Outputs: 
    
        Return a fasta file in the output directory.

    """
    # open the file to write the data in
    with open(out_file, 'w') as fout:
        # iterates to the dictionary items
        for name, seq in sequence_dict.items():
            # write the header
            fout.write(f'>{name}\n')
            # them write a line with 80 character at time util finish
            # the sequence len
            for i in range(0, len(seq), wrap):
                fout.write(f'{seq[i:i + wrap]}\n')


def fasta_reader(filename):
    """
    Read a multi or single fasta file.
    
    Inputs:
    
        filename - string that represents a name of the file or a path to
                   the file.
    
    Outputs:
    
        A generator object containing a Seq and ID biopython objects.
    """
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            for record in FastaIterator(handle):
                yield str(record.id), str(record.seq)
    else:
        with open(filename) as handle:
            for record in FastaIterator(handle):
                yield str(record.id), str(record.seq)                
