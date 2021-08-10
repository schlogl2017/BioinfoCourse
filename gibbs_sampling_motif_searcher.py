#usr/bin/env python
import math
from collections import defaultdict, Counter
from random import randint
from pprint import pprint


def most_frequent_patterns(sequence, k, n=1):
    """Returns the n most frequent patterns of lenght k from a input
    sequence."""
    return Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)]).most_common(n)


def get_random_motifs(motifs, k):
    """
    Function to generate random motifs from a 
    matrix (list of lists) of motifs (substrings with
    length m).
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        k - a integer representing the length of the random motifs
            expected to be generated by the function.
    
    Outputs:
        rand_motifs - an array/list-like object with the random
                      motifs of length k generated by the function.
    """
    # gets the motif length and the number of motifs
    ml, n = len(motifs[0]), len(motifs)    
    # chacks if all motifs have the same length
    assert all(len(m) == ml for m in motifs), print('The motifs must have the same length')     
    # initialize the container for thr random motifs
    rand_motifs = []
    # iterates through indexes of the motifs array
    for i in range(n):
        # calculates the random index
        ii = randint(0, (ml - k))
        # get the motif i from the array
        # and make a random motif with the random
        # indexes
        mot = motifs[i][ii:ii+k]
        # add the new motif to the container
        rand_motifs.append(mot)
    return rand_motifs


def get_most_probable_Kmer_from_profile(motifs, k, profile):
    """
    Function the find the most probable kmer of length k from a matrix of
    motifs.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.    
        k - a integer representing the length of the random motifs
            expected to be generated by the function.
        profile - a matrix-like object that represents the frequency of each
                  base from the allowed alphabet that represent a motif.  
    
    Outputs:
        prob_kmer - a string that represents the most probable kmer given a
                    profile. 
    """
    # get motif length
    ml = len(motifs[0])
    # chacks if all motifs have the same length
    assert all(len(m) == ml for m in motifs), print('The motifs must have the same length')     
    # get the number of motifs
    n = len(motifs)
    # initialize the counter
    # must be -1 because some kmer my have zero probability
    prob = -1
    # start with any k motifs
    prob_kmer = motifs[:k]
    # get the range of motif indexes
    for i in range(n - k + 1):
        # get the kmer to profile
        seq = motifs[i:i+k]
        # get the kmers probability
        prob_most = motif_probability(seq, profile)
        # if the kmer probability large than the
        # probability in the counter
        if prob_most > prob:
            # then it will be the new prob value
            prob = prob_most
            # and the kmer with the large probability
            # is the most probable one
            prob_kmer = seq
    return prob_kmer


def consensus_sequence(motifs, alphabet):
    """
    Function that finds the consensus sequence from a matrix of 
    motifs.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'         
    
    Outputs:
        consensus - a string that representes the most frequent bases
                    found in each motif position in the motif matrix.
    """
    # get the motif length and the number of motifs 
    k, n = len(motifs[0]), len(motifs)
    # chacks if all motifs have the same length
    assert all(len(m) == k for m in motifs), print('The motifs must have the same length')
    # count the bases at each position in the motifs
    counts = motifs_count(motifs, alphabet)
    # get the bases that belongs to the consensus
    # according to counts
    consensus = []
    # iterates through the motif length
    for i in range(k):
        # initialize the count to compare
        count = 0
        # the base container
        freq_base = ''
        # get the base from alphabet
        for base in alphabet:
            # checks if the base has a count bigger tha count
            if counts[base][i] > count:
                # if true counts get the base count
                count = counts[base][i]
                # then the current base is the most
                # frequent for tha particular postion in the motif
                freq_base = base
        # then the current base goes to consensus
        consensus.append(freq_base)
    return ''.join(consensus)



def consensus_sequence_with_pseudo(motifs, alphabet, pseudo=1):
    """
    Function that finds the consensus sequence from a matrix of 
    motifs.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'         
    
    Outputs:
        consensus - a string that representes the most frequent bases
                    found in each motif position in the motif matrix.
    """
    # get the motif length 
    k = len(motifs[0])
    # chacks if all motifs have the same length
    assert all(len(x) == k for x in motifs), print('The motifs must have the same length')
    # count the bases at each position in the motifs
    counts = motif_counts_pseudo(motifs, alphabet, pseudo)
    # get the bases that belongs to the consensus
    # according to counts    
    consensus = []
    # iterates through the motif length
    for i in range(k):
        # initialize the count to compare
        count = 0
        # the base container
        freq_base = ''
        # get the base from alphabet
        for base in alphabet:
            # checks if the base has a count bigger tha count
            if counts[base][i] > count:
                # if true counts get the base coun
                count = counts[base][i]
                # then the current base is the most
                # frequent for tha particular postion in the motif                
                freq_base = base
        # then the current base goes to consensus        
        consensus.append(freq_base)
    return ''.join(consensus)


def motifs_count(motifs, alphabet):
    """
    Returns the count of the nucleotides that appears in each
    columns of the motifs. Motifs are a matrix of strings(list of lists)
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'
    
    Outputs:
        count - a matrix-array-like that represent the count of the bases
                by position in a matrix of motifs.
    """
    # get the motif and the number of motifs 
    k, n = len(motifs[0]), len(motifs)
    # chacks if all motifs have the same length
    assert all(len(m) == k for m in motifs), print('The motifs must have the same length')         
    # initialize the counter
    count = defaultdict(int, [(base, [0] * k) for base in alphabet])
    # iterates through motifs
    for i in range(n):
        # iterates through motif length
        for j in range(k):
            # check the base at motif i in position j
            base = motifs[i][j]
            # then add one for that particular base
            count[base][j] += 1
    return count


def motifs_scores(motifs, alphabet):
    """
    Calculates the score of the motifs from a matrix of motifs.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'
    
    Outputs:
        score - a integer representing the score of the motif consensus
                generated by a motif matrix.    
    """
    # get the motif and the number of motifs 
    k, n = len(motifs[0]), len(motifs)      
    # chacks if all motifs have the same length
    assert all(len(m) == k for m in motifs), print('The motifs must have the same length')  
    # start the counter
    score = 0
    # get the consensus sequence from motifs
    consensus = consensus_sequence(motifs, alphabet)
    # iterates through motifs
    for i in range(n):
        # iterates through motif length
        for j in range(k):
            # compares bases at position j from motif and consensus
            # sequence
            if motifs[i][j] != consensus[j]:
                # if they are not the same
                # add one to counter
                score += 1
    return score


def motif_probability(motif, profile):
    """
    Calculates the probability that sequence fits the profile.

    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.    
        profile - a matrix-like object that represents the frequency of each
                  base from the allowed alphabet that represent a motif.
    
    Outputs:
        probability - a float the represent the probability that given motif
                      came from the given profile.
    """
    # start the counter
    probability = 1
    # iterate through the indexes and the bases of the motif
    for i, base in enumerate(motif):
        # multiply the counter with the 
        # base frequency in the profile
        probability *= profile[base][i]
    return probability


def entropy(sequence):
    """
    Calculates the entropy of a sequence.
    
    Inputs:
        sequence - a string that representes the sequence. 
    
    Outputs:
        a float number representing the sequence entropy.       
    """
    # get the bases counted and the length of the sequence
    p, lns = Counter(sequence), float(len(sequence))
    # calculates the sequence entropy
    return -sum( count/lns * math.log(count/lns, 2) for count in p.values())


def motif_counts_pseudo(motifs, alphabet, pseudo=1):
    """
    Function to count the bases in the motifs matrix
    adding a pseudo count to the cases where the count
    is zero.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'
    
    Outputs:
        counts - a dictionary-like object mapping the bases/string characters to their
                 calculated counts added a pseudo count.
    """
    # get the motif length and the total
    # number of motifs in  the array
    k, m = len(motifs[0]), len(motifs)
    assert all(len(m) == k for m in motifs), print('The motifs must have the same length')
    # initialize the counter with the pseudo count
    counts = defaultdict(int, [(base, [pseudo] * k) for base in alphabet])
    # iterates through the motifs indexes
    for i in range(m):
        # iterates through motif indexes
        for j in range(k):
            # gets the base in the motif i
            # postion j
            base = motifs[i][j]
            # adds one to the base in the positon j
            # in the counter
            counts[base][j] += 1
    return counts


def motif_profile_pseudo_counts(motifs, alphabet, pseudo=1):
    """
    Calculates the profile of a matrix of motifs.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'
        pseudo - a number representing the value of the pseudo count
                 to add to the base count in the motifs matrix.
                 Default is 1.
    Outputs:
        profile - a matrix-like object that represents the frequency of each
                  base from the allowed alphabet that represent a motif.    
    """
    # get the motif length
    k = len(motifs[0])
    # get the number of motifs
    n = len(motifs)    
    # chacking if the motifs have the same length
    assert all(len(m) == k for m in motifs), print('The motifs must have the same length')
    # count the base with pseudo counts
    counts = motif_counts_pseudo(motifs, alphabet, pseudo)
    # get the motif profiles
    profile = {base : [c / (n + 4) for c in values] for base, values in counts.items()}
    return profile
    
    
def gibbs_sampler_motif_search(sequences, alphabet, k, num_searches, pseudo):
    """
    Motif Discovery with Gibbs Sampling Algorithm that is a Probabilistic Optimization
    algorithm. 
    General definition: class of computational algorithms that rely on repeated
    random sampling to compute their results.
    Specific definition: randomized algorithm where the computational resources 
    used are bounded but the answer is not guaranteed to be correct 100% of the time.
    Function that return the best motifs given a matrix of sequences.
    
    Inputs:
        motifs - a matrix-like object (list of list) representing
                 a collection of substrings of length m.
        alphabet - a set of allowed character that represent the
                   alphabet that represents the bases in the motifs.
                   Ex: DNA = 'ACGT'
                       Protein = 'ACDEFGHIKLMNPQRSTVWY'
        k - a integer representing the length of the random motifs
            expected to be generated by the function.                       
        num_searches - integer representing how many times the search must be done
                       to select the best propable cabdidates motifs.
        pseudo - a real number representing the value of the pseudo count
                 to add to the base count in the motifs matrix.
                 Default is 1.    
    
    Outputs:
        best_motifs - a list/array-like object with the best putative cndidates motifs 
                      given a matrix of sequences.
    """
    # get the number of motifs
    n = len(sequences)
    # get the new random motifs
    motifs = get_random_motifs(sequences, k)
    # the current best motifs are the new random motifs
    best_motifs = motifs
    # iterates through the length of designed number o searchs
    # must be large to get stable results
    for j in range(1, num_searches):
        # get the random index
        i = randint(0, n - 1)
        # pops out from motifs a random selected motif
        motifs.pop(i)
        # get the profile of the motf matrix
        profile = motif_profile_pseudo_counts(motifs, alphabet, pseudo)
        # insert back the poped sequence and get the most probable kmer/motif
        # given the calculated profile
        motifs.insert(i, get_most_probable_Kmer_from_profile(sequences[i], k, profile))
        # scores the motifs and the best motifs
        if motifs_scores(motifs, alphabet) < motifs_scores(best_motifs, alphabet):
            # if the new random motifs have the lowest score
            # then they are the current best motifs
            best_motifs = motifs
    return best_motifs
    
    
    
if __name__ == "__main__":
      seq1="""atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccgacccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaagggggggatgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccgagctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggagatcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttataggtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcacccaaattcagtgtgggcgagcgcaacggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcataacttgagttaaaaaaaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgtattggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcataaaaaaaagggggggaccgaaagggaagctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga"""
      print(entropy(seq1))
      # 1.9747521324833461
      print(most_frequent_patterns(seq1, 10, n=1))
       
      m = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
      alphabet = 'ACGT'
      pprint(motifs_count(m, alphabet))
      # {'A': [1, 2, 1, 0, 0, 2],
      # 'C': [2, 1, 4, 2, 0, 0],
      # 'G': [1, 1, 0, 2, 1, 1],
      # 'T': [1, 1, 0, 1, 4, 2]})
        
      pprint(motif_counts_pseudo(m, alphabet, pseudo=1))
      # {'A': [2, 3, 2, 1, 1, 3], 
      # 'C': [3, 2, 5, 3, 1, 1], 
      # 'G': [2, 2, 1, 3, 2, 2], 
      # 'T': [2, 2, 1, 2, 5, 3]}
      
      pprint(motif_profile_pseudo_counts(m, alphabet, pseudo=1))
      #{'A': [0.2222222222222222, 
      #       0.3333333333333333, 
      #       0.2222222222222222, 
      #       0.1111111111111111,
      #       0.1111111111111111, 
      #       0.3333333333333333], 
      #'C': [0.3333333333333333, 
      #      0.2222222222222222, 
      #      0.5555555555555556, 
      #      0.3333333333333333, 
      #      0.1111111111111111, 
      #      0.1111111111111111], 
      # 'G': [0.2222222222222222, 
      #       0.2222222222222222, 
      #       0.1111111111111111, 
      #       0.3333333333333333, 
      #       0.2222222222222222, 
      #       0.2222222222222222], 
      # 'T': [0.2222222222222222, 
      #       0.2222222222222222, 
      #       0.1111111111111111, 
      #       0.2222222222222222, 
      #       0.5555555555555556, 
      #       0.3333333333333333]}
      
      pprint(consensus_sequence(m, alphabet)) 
      # CACCTA
        
      l2 = ['GGACTTCAGGCCCTATCGGA', 'GGTCCATCGACCCGCGGCCC', 'GGACGTAAGTCCCTAACGCG', 'CGATCATCGGCCAGGGCGCC', 
            'TGACCGACGTCCCCAGCCCC', 'GGACCTTCGGCCCCACCCAC', 'GGACTTCTGTCCCTAGCCCT', 'GGACTTTCGGCCCTGTCCGC', 
            'GGACTAACGGCCCTCAGGTG', 'GGACCGAAGTCCCCGGGCTC']
      
      pprint(consensus_sequence_with_pseudo(l2, alphabet, pseudo=1))
      # 'GGACCTACGGCCCTAGCCCC'
        
      pprint(motifs_count(m, alphabet))
      # {'A': [1, 2, 1, 0, 0, 2],
      # 'C': [2, 1, 4, 2, 0, 0],
      # 'G': [1, 1, 0, 2, 1, 1],
      # 'T': [1, 1, 0, 1, 4, 2]})
      
      pprint(motifs_scores(m, alphabet))
      # 14
      
      s = 'ACGGGGATTACC'  # 0.0008398080000000002
      profile =  { 'A': [ 0.2 ,0.2 ,0.0 ,0.0 ,0.0 ,0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
                   'C': [ 0.1 ,0.6 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.4 ,0.1, 0.2, 0.4, 0.6], 
                   'G': [ 0.0 ,0.0 ,1.0 ,1.0 ,0.9 ,0.9 ,0.1 ,0.0, 0.0 ,0.0 ,0.0, 0.0], 
                   'T': [ 0.7 ,0.2 ,0.0 ,0.0 ,0.1 ,0.1 ,0.0, 0.5 ,0.8, 0.7, 0.3, 0.4]}
      pprint(motif_probability(s, profile))
       
      dnas = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]
      pprint(get_random_motifs(dnas, 3))
      # something like TTA GTC ACG ACG GAG
       
      seqs = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
              'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
              'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
      k, N = 8, 100
      pprint(gibbs_sampler_motif_search(seqs, alphabet, k, N, 1))
      # AACGGCCA AAGTGCCA TAGTACCG AAGTTTCA ACGTGCAA
    
    
    
    
