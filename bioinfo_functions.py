#!usr/bin/env python
# -*- coding: utf-8 -*-


import re


def get_one_to_three_letter_aac(seq, acc_dict, sep='-'):
    """
    Function that gets a protein sequence and convert a one
    letter aminoacid to a three letter aminoacid representation.
    
    Inputs:
        seq - a string representing a protein sequence.
        aac_dict - a dictionary-like object mapping a 
                   one letter amino acid to a three letter
                   amino acid representation.
        sep - a string representing the separator.
              Default is '-'.
              
    Outputs:
        three_letter_seq - a string representing a protein sequence
                           with three letter amino acid representation.
    Ex:
    aa = {'A': 'Ala','R': 'Arg','N': 'Asn','D': 'Asp','C': 'Cys','Q': 'Gln',
          'E': 'Glu','G': 'Gly','H': 'His','I': 'Ile','L': 'Leu','K': 'Lys',
          'M': 'Met','F': 'Phe','P': 'Pro','S': 'Ser','T': 'Thr','W': 'Trp',
          'Y': 'Tyr','V': 'Val'}
    >> get_one_to_three_letter_aac('wvvafcc', aa, sep='-')
    'Trp-Val-Val-Ala-Phe-Cys-Cys'
    
    >>> get_one_to_three_letter_aac('', aa, sep='-')
    ''
    
    >>> get_one_to_three_letter_aac('xvvffaz', aa, sep='-')
    /tmp/ipykernel_3460/2339983473.py:32: UserWarning: X is not a valid amino acid in position 0
    warnings.warn(f"{acc} is not a valid amino acid in position {i}", UserWarning)
    /tmp/ipykernel_3460/2339983473.py:32: UserWarning: Z is not a valid amino acid in position 6
    warnings.warn(f"{acc} is not a valid amino acid in position {i}", UserWarning)
    '?-Val-Val-Phe-Phe-Ala-?'
    """
    three_letter_seq = []
    seq = seq.upper()
    for i, acc in enumerate(seq):
        # try to translante one letter to three letter
        try:
            three_letter_seq.append(acc_dict[acc])
        except:
            if acc not in acc_dict.keys():
                # if ambiguity in aac representation
                warnings.warn(f"{acc} is not a valid amino acid in position {i}", UserWarning)
                # add a question mark to sinalize the ambiguity
                three_letter_seq.append('?')
    if sep:
        return sep.join(three_letter_seq)
    return three_letter_seq
    

def get_codons(seq):
    seq = seq.upper()
    len_seq = len(seq)
    return [seq[i:i+3] for i in range(len_seq - 2, 3)]

























