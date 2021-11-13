#!/usr/bin/env python
# -*- coding: utf-8 -*-


def get_pur_pyr_counts(sequence):
    """
    Function to calculate the frequency of purines (A or G) and
    pyrimidines (T, C or U) in a sequence.

    Inputs:
        sequence - a stringt object representing a sequence (DNA,
        RNA).

    Outputs:
        pur - a integer number representing the total number of
              purines in a sequence.
        pyr - a integer number representing the total number of
              pyrimidines in a sequence.
    """
    pur = sequence.count('A') + sequence.count('G')
    pyr = sequence.count('T') + sequence.count('C') + sequence.count('U')
    return pur, pyr