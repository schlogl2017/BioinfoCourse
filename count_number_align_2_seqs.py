#!/usr/bin/env python

#I don't know what your m and n are, but keep in mind that this solution grows quickly. Also, I did some digging and scipy #approximates the calculation of the binomial coefficient, which is probably going to be an issue for large m and n.

import sys
import scipy.special

def alignments(m, n):
    if m < 0 or n < 0:
        raise ValueError('m and n should be non-negative')
    s = 0
    b = min(m, n) + 1
    for k in xrange(0, b):
        s += (2 ** k) * scipy.special.binom(m, k) * scipy.special.binom(n, k)
    return s

def main():
    m = 4
    n = 2
    a = alignments(m, n)
    assert (m == 4 and n == 2 and a == 41), "something went wrong!"
    sys.stdout.write("alignments(%d, %d) = %d\n" % (m, n, a))

if __name__ == "__main__":
    main()
