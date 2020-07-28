#!usr/bin/env python
# -*- coding: utf-8 -*-
import random
import gzip
import argparse
import numpy as np
import itertools
from collections import Counter
from alphabet import alphabet_dna


def is_header(line):
    return line[0] == '>'


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
    print(f'End reading fastq file: {filename}')


def nuc_count(sequence):
    sequence = sequence.upper()
    alt = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '-']
    return {base: sequence.count(base) for base in sequence if base not in alt}


def count_alternative_bases(sequence):
    alt = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '-']
    seq = sequence.upper()
    counter = {base: seq.count(base) for base in seq if base in alt}
    return counter


def nuc_frequency(sequence):
    alt_base = sum(count_alternative_bases(sequence).values())
    l_seq = len(sequence) - alt_base
    print(alt_base)
    counts = nuc_count(sequence)
    return {base: (value / l_seq) for base, value in counts.items()}


def count_kmers(sequence, k=1):
    return Counter([sequence[i:i+k].upper() for i in range(0, len(sequence) - k + 1)])


def kmer_frequency(sequence, k):
    l_seq = len(sequence)
    kmers = count_kmers(sequence, k)
    freq = {kmer: (count / l_seq) for kmer, count in kmers.items()}
    return freq


def create_random_sequence(alphabet, nuc_freqs, n):
    return ''.join(np.random.choice(alphabet, p=nuc_freqs) for _ in range(n))


def get_gc_content(sequence):
    """Returns the gc content of a genome."""
    len_seq = len(sequence) - sum(count_alternative_bases(sequence).values())
    sequence = sequence.upper()
    c = sequence.count('C')
    g = sequence.count('G')
    return round((c + g) / len_seq, 4)


def get_at_content(sequence):
    """Returns the at content of a genome."""
    return 1 - get_gc_content(sequence)


def get_at_gc_ratio(sequence):
    """Returns the at/gc ratio of a genome."""
    return get_at_content(sequence) / get_gc_content(sequence)


def get_strand_complement(sequence):
    """Returns the complement strand of the genome."""
    change = str.maketrans('ACGT', 'TGCA')
    return sequence.translate(change)


def get_reverse_complement(sequence):
    """Returns the reverse complement strand of the genome."""
    return get_strand_complement(sequence)[::-1]


def check_is_palindrome(sequence, kmer):
    """Returns True if the sequence is palindromic other wise
    False."""
    return sequence.find(kmer[::-1]) == 0


def get_palindromes(sequence, k):
    """Returns the count of all the palindromic
    substrings of a genome."""
    kmers = list(count_kmers(sequence, k).keys())
    rev_kmers = [get_reverse_complement(kmer) for kmer in kmers]
    palindromes = []
    for mer1, mer2 in zip(kmers, rev_kmers):
        if mer1 == mer2:
            palindromes.append((mer1, mer2))
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
    half = len(sequence) // 2
    genome = np.frombuffer(sequence.encode() + sequence.encode(), dtype='uint8')
    g = np.concatenate(([0], np.array(genome == ord('G'), dtype='uint8').cumsum()))
    c = np.concatenate(([0], np.array(genome == ord('C'), dtype='uint8').cumsum()))
    gc = g - c
    skew = gc[half:(half + len(sequence))] - gc[0:len(sequence)] + \
           gc[(len(sequence) - half):(2 * len(sequence) - half)] - gc[len(sequence):(2 * len(sequence))]
    return skew


def get_multinomial_model(sequence):
    """Returns a simple multinomial probabilistic model"""
    return nuc_frequency(sequence)


def get_markov_model(sequence):
    """Returns a simple Markov model"""
    alphabet = set(sequence)
    mkm = dict()
    for base in alphabet:
        intervals = sorted([0] + [random.random() for _ in range(3)] + [1])
        print(intervals)
        probs = [intervals[i+1] - intervals[i] for i in range(4)]
        mkm[base] = {base: prob for base, prob in zip(alphabet, probs)}
    return mkm


def main():
    parser = argparse.ArgumentParser('Get a basic genome statistics')
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('-p',
                        '--path',
                        type=str,
                        required=True,
                        dest='filename',
                        action='store',
                        help='The filename')
    parser.add_argument('-o',
                        '--outputfile',
                        dest='output',
                        help='Name of the output file')
    parser.add_argument('-w',
                        '--window',
                        type=int,
                        dest='window_size',
                        help='Length of the slide window')
    parser.add_argument('-s',
                        '--step',
                        type=int,
                        dest='step',
                        default=1,
                        help='Length of the step to be used')
    parser.add_argument('-k',
                        type=int,
                        dest='k',
                        default=1,
                        help='Length of the kmer to be searched')
    parser.add_argument('-m',
                        type=str,
                        dest='multinomial',
                        help='To return a multinomial probabilistic model')
    parser.add_argument('-mk',
                        type=str,
                        dest='markov',
                        help='To return a markov probabilistic model')
    parser.add_argument('-pal',
                        type=str,
                        dest='markov',
                        help='To return a palindromic analysis')
    parser.add_argument('-s',
                        '--skew',
                        type=str,
                        dest='skew',
                        help='To return an gc skew analysis')
    parser.add_argument('-r',
                        '--rand',
                        type=str,
                        dest='random',
                        help='To make a randomic sequence')
    parser.add_argument('-pr',
                        '--probabilities',
                        nargs=6,
                        metavar=('a', 'b', 'c', 'd'),
                        help="List of probabilities",
                        type=float,
                        default=None)
    parser.add_argument('-rv',
                        '--rev_comp',
                        help="To get the reverse complement of the sequence",
                        type=str)
    parser.add_argument('-gc',
                        type=str,
                        dest='gc_content',
                        help='To make a gc content analysis of a sequence')
    parser.add_argument('-at',
                        type=str,
                        dest='at_content',
                        help='To make a gc content analysis of a sequence')
    parser.add_argument('--at_ratio',
                        type=str,
                        dest='ratio',
                        help='To make a gc content analysis of a sequence')
    args = parser.parse_args()
    process_file(args)


def process_file(args):
    if args.window_size and args.step:
        for name, sequence in simple_fasta_reader(args.filename):
            window = get_chunks(sequence, args.window_size, args.step)
            gc_cont = [(get_gc_content(seq), end) for (seq, end) in window]
    elif args.k and args.window_size or args.step:
        for ids, sequence, quals in simple_fasta_reader(args.filename):
            window = get_chunks(sequence, args.window_size, args.step)
            kmers = [(count_kmers(seq), end) for (seq, end) in window]
            k_freq = kmer_frequency(sequence, args.k)
    elif args.k and args.pal:
        for ids, sequence, quals in simple_fasta_reader(args.filename):
            window = get_chunks(sequence, args.window_size, args.step)
            kmers = [(count_kmers(seq), end) for (seq, end) in window]
            pal = get_palindromes(sequence, args.k)
            k_freq = kmer_frequency(sequence, args.k)


if __name__ == '__main__':
    main()