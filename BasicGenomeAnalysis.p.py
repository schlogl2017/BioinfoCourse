#!usr/bin/env python
# -*- coding: utf-8 -*-
from collections import Counter, defaultdict


class GenomicSequence:
    """Represents a genome sequence"""

    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.length = len(sequence)
        self.bases = self.get_base_count
        self.base_composition = self.get_base_composition
        self.gc_content = self.get_gc_content()
        self.at_content = self.get_at_content()
        self.sequence_range = self.get_sequence_range
        # self.sequence_kmers = self.get_kmers()

    def __str__(self):
        return self.sequence

    def __getitem__(self, index):
        return self.sequence[index]

    def __index__(self):
        return 1

    def get_slice(self, start, end, pace=1):
        """Returns the substring with start position, end position
        and a strider."""
        return self.sequence[start:end:pace]

    def get_base_count(self):
        """Returns the count of the characters that compound the
        input string"""
        return Counter(self.sequence)

    def get_base_composition(self):
        """Returns the frequency of characters that compound the
        input string."""
        base_comp = {base: self.sequence.count(base)/self.length
                     for base in self.sequence}
        return base_comp

    def get_gc_content(self):
        """Returns the gc content of a genome."""
        c = self.sequence.count('C')
        g = self.sequence.count('G')
        return round((c + g) / self.length, 4)

    def get_at_content(self):
        """Returns the at content of a genome."""
        return 1 - self.get_gc_content()

    def get_at_gc_ratio(self):
        """Returns the at/gc ratio of a genome."""
        return self.get_at_content()/self.get_gc_content()

    def get_sequence_range(self, start, end, pace):
        """Returns the substring with start position, end position
                and a strider."""
        return self.sequence[start:end:pace]

    def get_kmers_counts(self, k=1):
        """Returns the count of all the contiguous and overlapping
        substrings of length K from a genome."""
        return Counter(self.sequence[i:i+k] for i in range(self.length - k + 1))

    def get_kmers_frequencies(self, k=1):
        """Returns the frequencies of all the contiguous and overlapping
                substrings of length K from a genome."""
        kmers = self.get_kmers_counts(k)
        freq = defaultdict(float)
        for mer, count in kmers.items():
            freq[mer] = count / sum(kmers.values())
        return freq

    def get_strand_complement(self):
        """Returns the complement strand of the genome."""
        change = str.maketrans('ACGT', 'TGCA')
        return self.sequence.translate(change)

    def get_reverse_complement(self):
        """Returns the reverse complement strand of the genome."""
        return self.get_strand_complement()[::-1]

    def check_is_palindrome(self, sequence):
        """Returns True if the sequence is palindromic other wise
        False."""
        return self.sequence.find(sequence[::-1]) == 0

    def get_palindromes(self, k):
        """Returns the count of all the palindromic
                substrings of a genome."""
        kmers = list(self.get_kmers_counts(k).keys())
        rev_kmers = [self.kmer.get_reverse_complement() for self.kmer in kmers]
        palindromes = []
        for mer1, mer2 in zip(kmers, rev_kmers):
            if mer1 == mer2:
                palindromes.append((mer1, mer2))
        return palindromes

    def get_palindromes_counts(self, k):
        """"Returns the count of all palindromes from a string."""
        palindromes = self.sequence.get_palindromes(k)
        return  Counter(palindromes)


# dna = GenomicSequence('aaaacgtagcc')
# print(dna)
# print(dna.get_base_count())
# print(dna.get_base_composition())
# print(dna.get_gc_content())
# print(dna.get_at_content())
# print(dna.get_slice(0, 6, 2))
# print(dna.get_sequence_range(0, 6, 2))
# print(dna.get_kmers_counts(3))
# print(dna.get_kmers_frequencies(3))
# print(sum(dna.get_kmers_frequencies(3).values()))
# print(dna.get_strand_complement())
# print(dna.get_reverse_complement())
# s = GenomicSequence('AGA')
# s2 = "AG"
# print(s.check_is_palindrome(s2))