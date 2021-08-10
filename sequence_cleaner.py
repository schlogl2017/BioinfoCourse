#!usr/bin/env python
# -*- coding: utf-8 -*-


def seq_cleaner(sequence, alphabet):
    """
    Clean up a sequence from not allowed characters.
    Input:
        sequence - sequence or a string
    Output:
        sequence - cleaned sequence or a string
    """
    seq = sequence.upper()
    sequence = [base for base in seq if base in alphabet]
    return ''.join(sequence)