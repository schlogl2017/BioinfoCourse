#!/usr/bin/env python
# coding: utf-8
import argparse
import gzip
import os
import time
import numpy as np
import itertools
from collections import Counter, defaultdict
import matplotlib.pyplot as plt

start_time = time.process_time()

parser = argparse.ArgumentParser(prog='Kmer analysis',
                                 usage='%(prog)s [options] path',
                                 description='Search for all k mers of length k from a genome.')
parser.add_argument('--path',
                    '-p',
                    metavar='path',
                    type=str,
                    required=True,
                    dest='path',
                    help='Path to the files')
parser.add_argument('--out', '-o',
                    type=str,
                    dest='output',
                    help='name output files.')
parser.add_argument('--counts',
                    '-c',
                    type=str,
                    store=True,
                    dest='count',
                    help='Count kmers')
parser.add_argument('--freq',
                    '-f',
                    type=str,
                    store='store_true',
                    dest='k',
                    help='Frequency of the kmers')
parser.add_argument('--pal',
                    '-p',
                    type=str,
                    store='store_true',
                    dest='k',
                    help='Check for palindromic sequences.')
parser.add_argument('--plot',
                    '-plt',
                    type=str,
                    store='store_true',
                    dest='k',
                    help='Check for palindromic sequences.')
args = parser.parse_args()


alphabet_dna = ['A', 'C', 'G', 'T']


def is_header(line):
    """Check if the line starts with '>'."""
    return line[0] == '>'


def parse_fasta_file(filename):
    """It reads and returns a name and the sequence from a file (compressed or not)."""
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        fasta_iter = (it[1] for it in itertools.groupby(f, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequences = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequences


def get_strand_complement(sequence):
    """Returns the complement strand of the genome."""
    seq = sequence.upper()
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """Returns the reverse complement strand of the genome."""
    seq = sequence.upper()
    return get_strand_complement(seq)[::-1]


def create_random_sequence(alphabet, nuc_freqs, n):
    return ''.join(np.random.choice(alphabet, p=nuc_freqs) for _ in range(n))


def kmer_count(sequence, alphabet, k=1):
    seq = sequence.upper()
    seq_len = len(seq)
    kmers = [seq[i:i + k] for i in range(0, seq_len - k + 1)]
    filterd_kmers = [kmer for kmer in kmers if all(base in set(alphabet) for base in kmer)]
    return Counter(filterd_kmers)


def print_kmers_stats(sequence, alphabet, start, k):
    seq = sequence.upper()
    print('  k     k-mers              4^k      N-k+1          missing   repeated')
    for k in range(start, k):
        kmers = kmer_count(seq, alphabet, k)
        print(
            f"{k} {len(kmers)} {4 ** k} {(len(seq) - 1) - k + 1} {4 ** k - len(kmers)} {(len(seq) - 1) - k + 1 - len(kmers)}")
    return


def kmers_frequency(kmers_dict):
    num_kmers = len(kmers_dict.keys())
    kmer_freq = defaultdict(float)
    for kmer, count in kmers_dict.items():
        kmer_freq[kmer] = (count / num_kmers)
    return kmer_freq


def kmer_positions(sequence, k):
    """Returns the position of all k-mers and it complement pair from sequence as a dictionary"""
    seq = sequence.upper()
    kmer_position = {}
    for i in range(1, len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmer_position[kmer] = kmer_position.get(kmer, []) + [i]
    # combine kmers with their reverse complements
    pair_position = {}
    for kmer, pos_lst in kmer_position.items():
        rev_mer = get_reverse_complement(kmer)
        if kmer < rev_mer:
            pair_position[kmer] = sorted(pos_lst + kmer_position.get(rev_mer, []))
        elif rev_mer < kmer:
            pair_position[rev_mer] = sorted(kmer_position.get(rev_mer, []) + pos_lst)
        else:
            pair_position[kmer] = pos_lst
    return pair_position


def check_is_palindrome(sequence, kmer):
    """Returns True if the sequence is palindromic other wise
    False."""
    seq = sequence.upper()
    kmer = kmer.upper()
    return seq.find(kmer[::-1]) == 0


def get_palindromes(sequence, alphabet, k):
    """Returns the count of all the palindromic
    substrings of a genome."""
    kmers = list(kmer_count(sequence, alphabet, k).keys())
    rev_kmers = [get_reverse_complement(kmer) for kmer in kmers]
    palindromes = set()
    for mer1, mer2 in zip(kmers, rev_kmers):
        if mer1 == mer2:
            palindromes.add((mer1, mer2))
    return palindromes


def get_chunks(sequence, window_size, step=1):
    """Returns a chunk of length of window_size and the end of the window size"""
    k = len(sequence)
    for i in range(0, k - window_size + 1, step):
        end = i + window_size
        chunk = sequence[i:i + window_size]
        assert len(chunk) == window_size
        yield chunk, end


def kmers_clumps(sequence, k, window, times):
    """ Find clumps of repeated k-mers in string. A clump occurs when t or more k-mers appear
        within a window of size w. A list of (k-mer, position, count) tuples is returned
        clumpList = kmers_clumps(seq, 9, 500, 6)
        print(len(clumpList))
        print([clumpList[i] for i in range(min(20,len(clumpList)))])
    """
    seq = sequence.upper()
    clumps = []
    kmers = kmer_positions(seq, k)
    for kmer, pos in kmers.items():
        for start in range(len(pos) - times):
            end = start + times - 1
            while (pos[end] - pos[start]) <= window - k:
                end += 1
                if end >= len(pos):
                    break
            if end - start >= times:
                clumps.append((kmer, pos[start], end - start))
    return clumps


def get_top_n_kmers(sequence, k, window, times):
    seq = sequence.upper()
    kmers = kmers_clumps(seq, k, window, times)
    # find kmers appearing in the most clumps
    mer_clump = {}
    for kmer, start, clumpSize in kmers:
        mer_clump[kmer] = mer_clump.get(kmer, 0) + 1
    top10 = [k for k in sorted(mer_clump, reverse=True, key=mer_clump.get)][:9]
    return top10


def plot_kmer_frequencies(sequence, start, kmers, top_10_kmers):
    plt.figure(num=None, figsize=(16, 6), dpi=100, facecolor='w', edgecolor='k')
    plt.plot([start, start], [0, 10], 'r--')
    for n, kmer in enumerate(top_10_kmers):
        positions = kmers[kmer]
        print(positions)
        plt.text(len(sequence), n + 0.4, kmer, fontsize=10)
        plt.plot(positions, [n + 0.5 for _ in range(len(positions))], 'o', markersize=4.0)
    plt.xlim((0, len(sequence)))
    plt.title('Most frequent k-mers')
    plt.xlabel('Sequence length')
    plt.ylabel('Number sequences')
    plt.tight_layout()

end = time.process_time() - start_time
print(f'This script take {end} seconds to finish')
print('Done')
