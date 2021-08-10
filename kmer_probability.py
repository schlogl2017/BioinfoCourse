#!usr/bin/env python
from collections import defaultdict
import toolz as tz



def string_probability(string, char_freqs):
    """
    Function to calculate the Bernoulli probability of a
    random string.
    
    Inputs:
        
        string -  a string object representing a word.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a float representing the probability of the input
               string. (IID - Independent and identically distributed
               of the charcaters).
    """
    # initialyze the counts
    prob = 1
    # iterates through the string charcaters
    for i in range(len(string)):
        # multiply the initial prob with the frequency
        # of the chars in the bacground char_freqs
        prob *= char_freqs[string[i]]
    return prob


def get_string_frequency(kmer_list, kmer_counts):
    """
    Function to calculate the frequency of a list of strings
    with length k.
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        kmer_counts - a dictionary-like object representing the
                     mapping of counts of the overlapping substring in a bigger string.
    
    Outputs:
    
        freqs - a float representing the frequencies of the input
               list of substring, found in a bigger strings
    """
    # initialize the countainer
    freqs = defaultdict(float)
    # get the toal number of substrings
    total = sum(kmer_counts.values())
    # iterates trhough the length of the 
    # substrings of interest of length k
    for kmer in kmer_list:
        # gets the frequency and adds to the countainer
        f = kmer_counts[kmer] / total
        freqs[kmer] = freqs.get(kmer, 0.0) + f
    return freqs


def kmer_list_probabilities(kmer_list, char_freqs):
    """
    Function to calculate the Bernoulli probability of a
    list of strings of length k (k-mers/n-grams).
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a dictionary-like object representing the probability of the input
               list of string. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    kmer_probs = defaultdict(float)
    for kmer in kmer_list:
        p = string_probability(kmer, char_freqs)
        kmer_probs[kmer] = kmer_probs.get(kmer, 0.0) + p
    return kmer_probs
    
    
def get_groups_probabilities(base_freqs, kmer_counts):
    """
    Function that groups the kmer with length betwen kmin (kmax - 2)
    and kmax by kmer lengths (uses toolz groupby). Then calculates
    the probability of each kmer groups based in the nucleotides
    sequence background.
    
    Inputs:
    
        kmer_counts - a dictionary-like object representing the mapping 
                      of counts of the overlapping substring in a bigger string.
        base_freqs - a dictionary-like object representing the mapping of the
                     nucleotides counted in a string to their caounts.
    
    Outputs:
    
        prob - a dictionary-like object mapping the kmers groups to their 
               probabilities. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    # initialize the dictionary of dictionaries
    probs = defaultdict(lambda: defaultdict(dict))
    # group the kmer by length (kmin-2 / kmax)
    kmer_groups = tz.groupby(len, kmer_counts)
    # get list of group keys
    group_keys = list(kmer_groups.keys())
    # iterates through the keys
    for k in group_keys:
        # iterates through the kmers groups
        for kmer in kmer_groups[k]:
            # calculates the probabilities for each kmer in each group
            p = string_probability(kmer, base_freqs)
            # add the key groups to the out dictionary
            # and kmer keys and their probabilities to the inner dictionary
            probs[k][kmer] = p
    return probs    
    
    
def kmer_list_probabilities(kmer_list, char_freqs):
    """
    Function to calculate the Bernoulli probability of a
    list of strings of length k (k-mers/n-grams).
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a dictionary-like object representing the probability of the input
               list of string. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    kmer_probs = defaultdict(float)
    for kmer in kmer_list:
        p = string_probability(kmer, char_freqs)
        kmer_probs[kmer] = kmer_probs.get(kmer, 0.0) + p
    return kmer_probs    
    
    
def get_count_frequencies(counts):
    """Receive a dictionary as key:counts."""
    total = sum(counts.values())
    return {key: (cnt/total) for key, cnt in counts.items()}   
    
                
def probability(length, mismatches, p_match) :
    """ this function is used to work out the probability of a palindrome of
    length - 'length', with x amount of mismatches, based on the probability of
    getting a match (p_match)

    the combination funtctio (comb) is used instead of scipy to work out the number
    of possible ways of getting that mismatch and length combination to work out
    the probability of that occuring at random

    probability is defined by:
    p = combinatons of mismatche (N choose K) * (probability of a match based on the
    codon bias ** number of matches) * (probability of a mismatch ** number of
    mismatches)

    (** to the power of)
    """
    half_len = length //2
    matches =  half_len - mismatches
    comb_mismatch = comb(half_len, mismatches)
    p_mis = 1 - p_match
    return comb_mismatch * ((p_match ** matches) * (p_mis ** mismatches))       
    
    
def sequence_probability(seq, base_composition):
    """
    The function uses the back ground composition for the sequence 
    in question to calculate the probability for a match (P(match)) and mismatch (P(mismatch)).
    A match is defined when a base from the 3’ end is equal to its reverse complement 
    in its equivalent position from the 5’ end.
    
    Inputs:
    
        seq - a string representin a sequence of some length. The given sequence is
              treated as a pallindrome and the function looks to find if its 
              opposite palindromic pair (the other half) is a match or mismatch, 
              once it looks for half length of the sequence.
        base_composition - a dictionary-like mapping the bases that are counted in the
                           sequence to their frequencies.
    
    Outputs:
        
        probability - probability of the sequence is occuring at random. 
        mismatches - number of bases mismatches beteewn the half length sequence.
    """
    # only need to use half of the length of the sequence
    half_len = len(seq)//2
    # get the mismatches counted in the sequence
    mismatches = count_sequence_mismatches(seq)
    # get the number of matches in the sequence
    matches = half_len - mismatches
    # the number of ways of getting that num of mismatches in a given sequence
    # comb is the scipy.stats function
    comb_mismatch = comb(half_len, mismatches)
    # check for the matching bases
    matching_bases = ['AT','TA', 'GC', 'CG']
    # get the matches probability
    p_match = sum([base_composition[b[0]] * base_composition[b[1]] for b in matching_bases])
    # get the mismatches probabilities
    p_mis = 1 - p_match
    # get the probability of a given sequence occuring
    # randomly in a sequence.
    probability = comb_mismatch * ((p_match ** matches) * (p_mis ** mismatches))
    # returns the 
    return probability, mismatches     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
