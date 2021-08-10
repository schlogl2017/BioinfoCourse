#!usr/bin/env python
from collections import defaultdict
from scipy.stats import binom

# The probability of the word can be estimated as the product of probability of its letters.

def word_probability(word, letters_frequencies):
    """
    Function to calculate the probability of a word.
    The probability is the multiplication of the frequency or
    occurencies of the letters in a text or big string.
    
    Inputs:
        word - a substring composed with letter with any length.
        letters_frequencies - a dictionary-like object mapping the
                              letter to their frequency in a context
                              of a text or a large string.
    
    Outputs:
        prob - a float number representing the calculated
        probability of the word in a text or a large string.
    """
    # initialize the counter
    prob = 1
    # iterates through each letter in the word
    for letter in word:
        # multiply each letter frequency to obtain
        # the word probability
        prob *= letters_frequencies[letter]
    return prob

# N = number of possible position in a text or large string
# that the word can be found
N = sys.argv[]
k = len(word)
postions = N - k + 1

# The expected number of occurrences is the product of the number of # possible positions by the probability of occurrence at a given 
# position.
prob = word_probability(word, letters_frequencies)
expected = postions * prob

# The probability of occurrences can be obtained with the binomial 
# distribution via a binomial function.
# Function receive two parameters:
# n = number of trials 
# p =  is the probability of a single success
occurencies = binom.pmf(n = positions, p = prob)
