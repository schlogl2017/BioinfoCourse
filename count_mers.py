#!usr/bin/env python
# -*- coding: utf-8 -*-


import gzip
from collections import Counter
import itertools


def kmer_counter(sequence, k):
    return Counter(sequence[i:i+k] for i in range(0, len(sequence)-k+1))


def is_header(line):
    return line[0] == '>'


def count_fasta_files(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        return sum(g for g, _ in itertools.groupby(f, key=is_header))


def str_punctuation_strip(word):
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for _ in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


def simple_fasta_reader(filename):
    print('Starting reading the fasta file' + '\n')
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as fh:
        fasta_iter = (it[1] for it in itertools.groupby(fh, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequence = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequence


# print(kmer_counter("+TAGACAT",3))
# print(kmer_counter("+missmississippi",3))
for name, seq in simple_fasta_reader('VibrioCholerae.fa'):
    sequence = seq

print(' k       kmer            4^k         n-k+1           missing         repeated')
for k in range(3, 25):
    kmers = kmer_counter(sequence, k)
    print(f" {k}        {len(kmers)}        {4 ** k}        {(len(seq) - 1) - k + 1}        {4 ** k - len(kmers)}       {(len(seq) - 1) - k + 1 - len(kmers)}")