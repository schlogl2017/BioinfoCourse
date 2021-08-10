#!/usr/bin/env python
# coding: utf-8
import argparse
import gzip
import os
import time
import numpy as np
import itertools
from collections import Counter, defaultdict
from more_itertools import windowed
import matplotlib.pyplot as plt

start_time = time.process_time()

parser = argparse.ArgumentParser(description='Kmers from proteins.')
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
                    dest='mol_type',
                    help='Type mol DNA/RNA')
parser.add_argument('--length',
                    '-l',
                    type=int,
                    dest='k',
                    help='Length of the kmers')
args = parser.parse_args()


alphabet_dna = ['A', 'C', 'G', 'T']
alternative = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']


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


def nuc_count(sequence, alphabet):
    sequence = sequence.upper()
    return Counter(base for base in sequence if base in alphabet)


def alternative_bases_counter(sequence, alternative):
    seq = sequence.upper()
    return Counter(base for base in seq if base in alternative)


def nuc_frequency(sequence, counts, alternative):
    ns = alternative_bases_counter(sequence, alternative)
    l_seq = len(sequence) - sum(ns.values())
    return {base: (value / l_seq) for base, value in counts.items()}


def get_gc_content(sequence):
    """Returns the gc content of a genome."""
    len_seq = len(sequence) - sum(alternative_bases_counter(sequence).values())
    sequence = sequence.upper()
    c = sequence.count('C')
    g = sequence.count('G')
    return round((c + g) / len_seq, 4)


def get_at_content(gc):
    """Returns the at content of a genome."""
    return 1 - gc


def get_at_gc_ratio(at, gc):
    """Returns the at/gc ratio of a genome."""
    return at / gc


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


def kmers(sequence, alphabet, k):
    """Returns a generator of all mers(substring) of length k with overlap window
    from a string."""
    mers =  (''.join(c) for c in windowed(k, sequence))
    return [mer for mer in mers if all(base in set(alphabet) for base in mer)]


def kmers_counts(kmer_list):
    return Counter(kmer_list)


def check_is_palindrome(mer1, mer2):
    """Returns True if the sequence is palindromic other wise
    False."""
    return mer1.find(mer2[::-1]) == 0


def get_palindromes(kmer_list):
    """Returns the count of all the palindromic
    substrings of a genome."""
    rev_kmers = [get_reverse_complement(kmer) for kmer in kmer_list]
    palindromes = set()
    for mer1, mer2 in zip(kmer_list, rev_kmers):
        if check_is_palindrome(mer1, mer2):
            palindromes.add(mer1)
    return palindromes


def get_chunks(sequence, window_size, step=1):
    """Returns a chunk of length of window_size and the end of the window size"""
    k = len(sequence)
    for i in range(0, k - window_size + 1, step):
        end = i + window_size
        chunk = sequence[i:i + window_size]
        assert len(chunk) == window_size
        yield chunk, end


def gc_skew(sequence):
    """Finds the genome GC skew."""
    seq = sequence.upper()
    half = len(sequence) // 2
    genome = np.frombuffer(seq.encode() + seq.encode(), dtype='uint8')
    g = np.concatenate(([0], np.array(genome == ord('G'), dtype='uint8').cumsum()))
    c = np.concatenate(([0], np.array(genome == ord('C'), dtype='uint8').cumsum()))
    gc = g - c
    skew = gc[half:(half + len(sequence))] \
           - gc[0:len(sequence)] \
           + gc[(len(sequence) - half):(2 * len(sequence) - half)] \
           - gc[len(sequence):(2 * len(sequence))]
    return skew


def get_bases_stats(sequence, alphabet, start):
    """Returns the base statistics from different strands of a bacterial chromossome.
    It shows the difference in C e G composition in the lagg and lead strand"""
    seq = sequence.upper()
    seq_len = len(seq)
    half_seq = seq_len // 2
    ter = start + half_seq
    # as a circular genome
    if ter > seq_len:
        ter = ter - seq_len + 1
    counts = defaultdict(int)
    for base in alphabet:
        total = seq.count(base)
        if ter > start:  # start ---> ter
            f_count = seq[start:ter].count(base)
            r_count = total - f_count
        else:  # ter ---> start
            r_count = seq[ter:start].count(base)
            f_count = total - r_count
        counts[base] = (total, f_count, r_count)
    return counts


def print_strand_stats(sequence, alphabet, start):
    strands = get_bases_stats(sequence, alphabet, start)
    for base in alphabet_dna:
        total, f_count, r_count = strands[base]
        print(f'{base}:\t{total}\t{f_count}\t{r_count}\t{f_count - r_count}')
    return


def get_gc_strand_difference(sequence, start):
    seq = sequence.upper()
    seq_len = len(seq)
    half = seq_len // 2
    ter = start + half
    g = None
    c = None
    if ter > seq_len:
        ter = ter - seq_len + 1
    elif ter > start:
        g = 2 * seq[start:ter].count('G') - seq.count('G')
        c = 2 * seq[start:ter].count('C') - seq.count('C')
    else:
        g = seq.count('G') - 2 * seq[start:ter].count('G')
        c = seq.count('C') - 2 * seq[start:ter].count('C')
    diff = g - c
    return diff


def gc_skew(genome, n):
    x = []
    y = []
    for i in range(1, len(genome), n):
        x.append(i)
        diff = get_gc_strand_difference(genome, i)
        y.append(diff)
    return x, y


def get_sequence_skew(sequence):
    """Returns the difference between the total number of 
    occurrences of G and the total number of occurrences of C in 
    the first i elements of the sequence. """
    skew = [0]
    for idx, element in enumerate(sequence):
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            skew.append(skew[idx])
    return skew


def plot_gc_skew(x, y, start, offset):
    """Returns a plot of the genome gc skew
    input: GC_skew(sequence, n)"""
    plt.figure(num=None, figsize=(24, 7), dpi=100)
    yargmax = y.index(max(y))
    plt.axvline(start + offset, color="r", linestyle='--')
    plt.axvline(x[yargmax], color="g", linestyle='--')
    plt.plot(x, y)
    plt.savefig('Skew_genome.pdf', format='pdf', dpi=1200)
    plt.show()


def plot_base_counts_genomes(name, base_counts, length):
    base_markers = {'A': 'b-',
                    'C': 'r-',
                    'G': 'g-',
                    'T': 'y-',
                    'N': 'k-'}
    for base in ['A', 'C', 'G', 'T', 'N']:
        plt.plot(range(length), base_counts[base], base_markers[base], label=base)
    plt.xlabel('Position')
    plt.ylabel('Counts')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{name}_base_counts.png')
