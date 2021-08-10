#!usr/bin/env python
import numpy as np
from scipy.stats import poisson

def poisson_distribution(k, lambd):
    """
    Calculates the probability of observe k events
    with expectation mean lambd.
    
    Inputs:
        k - integer representing the k number of sucess.
        lambd - integer representing the expected mean
                for a sucess in a period of time.
    
    Outputs:
        A float number representing the probability of
        observe the occurence of k events.
    """
    return (lambd ** k * np.exp(-lambd)) / np.math.factorial(k)

def scipy_poisson(k, l):
    return poisson.pmf(k, l)
