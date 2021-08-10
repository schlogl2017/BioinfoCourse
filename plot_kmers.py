#usr/bin/env python
import matplotlib.pyplot as plt

def kmer_histogram(kmer_counts, k, max_counts):
    h = [0 for i in range(max_counts)]
    for kmer in kmer_counts:
        count = kmer_counts[kmer]
        if count < max_counts:
            h[count] += 1
    fig = plt.figure(figsize=(10, 6))
    plt.plot([i for i in range(max_counts)], h)
    plt.show();
