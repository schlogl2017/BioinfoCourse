# tatistical Analysis of Complete Bacterial Genomes: Avoidance
# of Palindromes and Restriction-Modification Systems

def contrast(kmer_obs, kmer_expec, len_seq, sigma):
    for kmer in kmer_obs:
        contrast = kmer_obs[kmer] - kmer_expec[kmer] / len_seq ** 1/2 * sigma
    return contrast

def variance(kmer_expec, kmer_count, len_seq, X):
    for kmer_expec:
    prefix = kmer[:-1]
    suffix = kmer[1:]
    middle = kmer[1:-1]
    var = (kmer_expec[kmer] / len_seq) * (1- (prefix/middle)) * (1-(suffix/middle))
    return var


def correlation(a, b, N, Da, Db):
    """
    a e b = contrasts from genome A and B.
    N = total number of m-letter words
    Da e Db = variations
    """
    cor = sum(a*b/(N-1) -sum(a)*sum(b)/(N-1)*N / Da ** 1/2 * Db ** 1/2
    return cor


def variations(N, a, b):
   Da = sum(a **b 2) / (N - 1) - (sum(a)) ** 2 / (N - 1) ** N
   Db = sum(b ** 2) / (N - 1) - (sum()) ** 2 / (N - 1) ** N
   return Da, Db  
   
   
def probability_not_find_kmer(len_seq, k):
    """
    Calculates the probability of not find the kmer in 
    any position in the sequence.
    
    Inputs:
        len_seq - integer representing the sequence length.
        k - integer representing the kmer length.
    
    Outputs:
        a float representing the probability of not find the kmer
        of length k in the sequence of length seq_length.
    """
    pos = (len_seq - k + 1)
    kmer_number = 1 / (4**k)
    return (1 - kmer_number) ** pos


def motif_probability(prob_not)
    """
    
    """
    return 1 - prob_not


def pmf(values):
    """
    values -- list of values
    Return probability mass function (pmf)
    dict key=value, value=prob
    """
    count = {}
    for v in values:
        count[v] = count.get(v, 0) + 1
    N = float(len(values))
    for v in count:
        count[v] /= N
    return count
