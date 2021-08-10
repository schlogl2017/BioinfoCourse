
from codon_translation import codon_map


def split_codons(orf):
    codon_list = []
    for i in range(0, len(orf)-2, 3):
        codon_list.append(orf[i:i+3])
    return codon_list


def count_codons(orf):
    orf = orf.upper()
    codon_list = split_codons(orf)
    counts = {}
    for codon in codon_map:
        counts[codon] = 0
    for codon in codon_list:
        counts[codon] += 1
    return counts


def group_counts_by_aa(codon_counts):
    counts_by_aa = {}
    for codon, acid in codon_map.items():
        if acid not in counts_by_aa:
            counts_by_aa[acid] = {}
        counts_by_aa[acid][codon] = codon_counts[codon]
    return counts_by_aa


def normalize_counts(grouped_counts):
    for aa in grouped_counts:
        total_count = sum(grouped_counts[aa].values())
        for codon in grouped_counts[aa]:
            grouped_counts[aa][codon] =  grouped_counts[aa][codon] / float(total_count)


def normalize_counts(grouped_counts):
    grouped_freqs = {}
    for aa in grouped_counts:
        grouped_freqs[aa] = {}
        total_count = sum(grouped_counts[aa].values())
        for codon in grouped_counts[aa]:
            grouped_freqs[aa][codon] = grouped_counts[aa][codon] / float(total_count)
    return grouped_freqs


def find_codon_bias(orf):
    codon_counts = count_codons(orf)
    stats_per_aa = group_counts_by_aa(codon_counts)
    return normalize_counts(stats_per_aa)
    

def print_codon_bias(orf):
    d = find_codon_bias(orf)
    for aa in d:
        print "%s:" % (aa)
        for codon in d[aa]:
            print "\t%s:\t%d%%" % (codon, int(d[aa][codon]*100))
            
if __name__ == "__main__":

    f = open('sample_orfs.txt', 'r')
    orf_list = list()
    for line in f:
        seq = line.strip()
        orf_list.append(seq)
    f.close()

    print_codon_bias(orf_list[0])

