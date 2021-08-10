#!usr/bin/env python


def get_all_mers_min_max(alphabet, kmin, kmax):
    """
    
    """
    # initialize the countainer
    all_kmers = []
    # iterates through the range kmin to kmax+1
    for n in range(kmin, kmax+1):
        # creates all kmers with length ranging from kmin to kmax
        kmers =  [''.join(letters) for letters in itertools.product(alphabet, repeat=n)]
        # append the list of kmers in all_kmers
        all_kmers.append(kmers)
    #list(itertools.chain(*all_kmers))
    return list(itertools.chain.from_iterable(all_kmers))


def get_all_mers(alphabet, kmin, kmax):
    all_kmers = []
    for n in range(kmin, kmax+1):
        kmer = [''.join(letters) for letters in itertools.product(alphabet, repeat=n)]
        all_kmers.append(kmer)
    for kmer in itertools.chain.from_iterable(all_kmers):
        yield kmer

def substitute(position, letter, string):
    return_value = ""
    if (position > 0):
        return_value = return_value + string[0:position]
    return_value = return_value + letter
    if (position < (len(string) - 1)):
        return_value = return_value + string[position+1:]                   
    return(return_value)

def make_upto_kmer_list(k_values, alphabet):

    # Compute the k-mer for each value of k.
    return_value = []
    for k in k_values:
        return_value.extend(make_kmer_list(k, alphabet))
    return(return_value)

def get_expected_kmer_zom(kmer, base_counts, len_sequence):
    n = len_sequence - len(kmer) + 1
    bases_kmer = list(kmer)
    expected = 1
    for base in bases_kmer:
        expected *= base_counts[base]
    return int(expected * n)
    
def expected_kmer_by_zom(kmer, base_freqs, len_seq):
    """
    Calculates the expected number of a substring of length k (k-mer)
    based in a zero order Markov model (zom).
    
    Inputs:
    
        kmer - a substring representing the word/k-mer (a string of length k).
        base_freqs - a dictionary-like object mapping the frequency of the
                     nucleotides/basesand their counted values normalized by
                     the length of the sequence/genome.
        len_seq - a integer representing the length of the sequence/genome where
                  the kmers were counted.
    
    Outputs:
        expected - a integer representing the kmer/substring of length k in the
                sequence/genome of length N (len_seq - len_kmer + 1).
    N = len(seq) - len(kmer) + 1
    E(w) = N*(nuc1*nuc2*nuc3)
    """
    # make a list of letter/kmer
    k_l = list(kmer)
    # number position in the genome where the
    # kmer was counted
    n = len_seq - len(kmer) + 1
    # counter to receive the bases values
    # to multiplicate
    cnt = 1
    # iterates in the kmer list
    # and gets the bases
    # to recover the frequencies from the dcitionary
    for base in k_l:
        # multiply all the bases values that
        # are found in the kmer
        cnt *= base_freqs[base]
    # return the expected value of the kmer
    return int(cnt * n)    
    
def count_all_bases(sequence):
    """
    Calculates the nucleotides/base composition of a DNA sequence.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.
    
    Outputs:
    
        all_bases - a dictionary-like object tha represent the count
                    of the nucleotides/bases that compound the sequence.
                    nucleotides are keys and the fruencies (as itegers) as values.
    
    """
    # create a set of bases
    bases = set(sequence)
    all_bases = defaultdict(int)
    # iterates in the base set
    for base in bases:
        # count the bases in the sequence
        all_bases[base] = sequence.count(base)
    return all_bases    
    
def expected_zero(kmer, len_seq, nuc_freqs):
    """
    Calculates the expected number of a substring of length k (k-mer)
    based in a zero order Markov model (zom).
    
    Inputs:
    
        kmer - a substring representing the word/k-mer (a string of length k).
        base_freqs - a dictionary-like object mapping the frequency of the
                     nucleotides/basesand their counted values normalized by
                     the length of the sequence/genome.
        len_seq - a integer representing the length of the sequence/genome where
                  the kmers were counted.
    
    Outputs:
        expected - a integer representing the kmer/substring of length k in the
                sequence/genome of length N (len_seq - len_kmer + 1).
    
    
    """
    # number positions where the kmer could be counted
    N = len_seq - len(kmer) + 1
    # get the base counts in the kmer sequence
    kbc = base_stats(kmer, 'ACGT', True, True)
    # the counter 
    expec = 1
    # iterates through the base kmer list
    # and the dictionary of nucleotides frequencies
    for bk, bs in zip(kbc, nuc_freqs):
        # multiply the nucleotides frequencies 
        # in the power of the base counts in the kmer
        expec *= nuc_freqs[bs]**kbc[bk]
    # returns the expected values fro the current kmer
    # in the genome
    return int(expec * N)    
    
def string_probability(string, char_freqs):
    """
    Function to calculate the Bernoulli probability of a
    random string.
    
    Inputs:
        
        string -  a string object representing a word.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a float representing the probability of the input
               string. (IID - Independent and identically distributed
               of the charcaters).
    """
    # initialyze the counts
    prob = 1
    # iterates through the string charcaters
    for i in range(len(string)):
        # multiply the initial prob with the frequency
        # of the chars in the bacground char_freqs
        prob *= char_freqs[string[i]]
    return prob    
    
def get_string_frequency(kmer_list, kmer_counts):
    """
    Function to calculate the frequency of a list of strings
    with length k.
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        kmer_counts - a dictionary-like object representing the
                     mapping of counts of the overlapping substring in a bigger string.
    
    Outputs:
    
        freqs - a float representing the frequencies of the input
               list of substring, found in a bigger strings
    """
    # initialize the countainer
    freqs = defaultdict(float)
    # get the toal number of substrings
    total = sum(kmer_counts.values())
    # iterates trhough the length of the 
    # substrings of interest of length k
    for kmer in kmer_list:
        # gets the frequency and adds to the countainer
        f = kmer_counts[kmer] / total
        freqs[kmer] = freqs.get(kmer, 0.0) + f
    return freqs    
    
def get_groups_probabilities(base_freqs, kmer_counts):
    """
    Function that groups the kmer with length betwen kmin (kmax - 2)
    and kmax by kmer lengths (uses toolz groupby). Then calculates
    the probability of each kmer groups based in the nucleotides
    sequence background.
    
    Inputs:
    
        kmer_counts - a dictionary-like object representing the mapping 
                      of counts of the overlapping substring in a bigger string.
        base_freqs - a dictionary-like object representing the mapping of the
                     nucleotides counted in a string to their caounts.
    
    Outputs:
    
        prob - a dictionary-like object mapping the kmers groups to their 
               probabilities. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    # initialize the dictionary of dictionaries
    probs = defaultdict(lambda: defaultdict(dict))
    # group the kmer by length (kmin-2 / kmax)
    kmer_groups = tz.groupby(len, kmer_counts)
    # get list of group keys
    group_keys = list(kmer_groups.keys())
    # iterates through the keys
    for k in group_keys:
        # iterates through the kmers groups
        for kmer in kmer_groups[k]:
            # calculates the probabilities for each kmer in each group
            p = string_probability(kmer, base_freqs)
            # add the key groups to the out dictionary
            # and kmer keys and their probabilities to the inner dictionary
            probs[k][kmer] = p
    return probs    
    
def from_kmer_list_to_probabilities(kmer_list, char_freqs):
    """
    Function to calculate the Bernoulli probability of a
    list of strings of length k (k-mers/n-grams).
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a dictionary-like object representing the probability of the input
               list of string. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    kmer_probs = defaultdict(float)
    for kmer in kmer_list:
        p = string_probability(kmer, char_freqs)
        kmer_probs[kmer] = kmer_probs.get(kmer, 0.0) + p
    return kmer_probs    
    
def expected_kmer_zero_markov(len_sequence, kmer_list, prob_dict):
    """
    Calculates the expected number of a kmer according to a 
    IID Bernoulli model(Zero Order Markov model). But here
    we calculate we the probabilities from a real world distribution
    what gives diferent expected values than the IID (A=T=C=G=0.25).
    
    Inputs:
    
        len_sequence - length of the sequence where the kmers were counted.
        kmer_list - list-like object representing all kmer counted in a sequence.
        prob_dict - dictionary-like object mapping the kmers to their calculated
                    probabilities.
    
    Outputs:
    
        kmer_expected - dictionary-like object mapping the kmers and their expected
                        distribution in a sequence of 'len_sequence' length.
    
    By definition the expected number of times we see an 'A' in our 'n' bp sequence in
    a random model (Bernoulli) is:
        EN = EX 1 + EX 2 + · · · + EX n = nEX 1 = 'npA'.
        
    >>> expected_kmer_zero_markov(len_sequence, ['GTTAAC'], p={A=T=.25; C=G=0.25})
    {'GTTAAC': 893}
    
    Calcuated based in a real counts in H.Influenza genome
    {'A': 0.3101722654062103, 'C': 0.19164984084916872,  'G': 0.1898531440827311,
    'T': 0.3083247496618899}
    >>> expected_kmer_zero_markov(len_sequence, ['GTTAAC'], p=0.0003327736695673287))
    {'GTTAAC': 1217}
    """
    # initialize the countainer for the probability data
    kmer_expected = defaultdict(int)
    # all kmes must have the same length
    k = len(kmer_list[0])
    # total length of the sequence
    N = len_sequence - k + 1
    # iterates over all kmers
    # get the probrabilities and multiply by the total
    # length of the sequence
    for kmer in kmer_list:
        exp = int(N * prob_dict[kmer])
        kmer_expected[kmer] = kmer_expected.get(kmer, 
                                                0) + exp
    return kmer_expected    
    
def get_expected_higher_markov(kmer_list, kmer_counts):
    """
    Calculates the expected value for a list of kmers and their counts.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
    
    Output:
    
        expected - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
    
    The expected values are calculated as:
    'Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]'
    """
    expected = defaultdict(int)
    for kmer in kmer_list:
        suf, pref, mid = kmer_counts[kmer[1:]], kmer_counts[kmer[:-1]], kmer_counts[kmer[1:-1]]
        if mid == 0:
            expected[kmer] = expected.get(kmer, 
                                          0)
        else:
            expected[kmer] = expected.get(kmer, 
                                          0) + int((pref * suf) / mid)
    return expected    
    
def get_expected_higher_markov(kmer_list, kmer_counts):
    """
    Calculates the expected value for a list of kmers and their counts.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
    
    Output:
    
        expected - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
    
    The expected values are calculated as:
    'Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]'
    """
    expected = defaultdict(int)
    for kmer in kmer_list:
        suf, pref, mid = kmer_counts[kmer[1:]], kmer_counts[kmer[:-1]], kmer_counts[kmer[1:-1]]
        if mid == 0:
            expected[kmer] = expected.get(kmer, 
                                          0)
        else:
            expected[kmer] = expected.get(kmer, 
                                          0) + int((pref * suf) / mid)
    return expected    
    
def get_variance(kmer_list, len_seq, kmer_expected):
    """
    Calculates the variance from a list of strings of length k.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_expectd - a dictionary-like object mapping kmer to their calculated
                       expected values.
    
    Outputs:
    
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
                   
    Because the model for the count is the sum of N almost independent observations, 
    each with probability P(W), it can be well modeled as a binomial distribution, 
    with varianceThe variance is calculated as:
    E(C(W)) * (1 - E(C(W))/N)    
    """
    k = len(kmer_list[0])
    N = len_seq - k + 1
    variance = defaultdict(float)
    for kmer in kmer_list:
        ex_val = kmer_expected[kmer]
        if ex_val == 0:
            variance[kmer] = variance.get(kmer, 0.0)
        else:
            var = ex_val * (1 - ex_val / N)
            variance[kmer] = variance.get(kmer, 0.0) + var
    return variance    
    
def get_standard_deviation(variance):
    """
    Calaculates the standard deviation from the kmers expected values.
    
    Inputs:
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
    
    Outputs:
        
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    The variance is calculated as:
    sigma(W) = sqrt(Expected) * (1 - Expected/len(seq) -k + 1))
    """
    # initialize the container
    std = defaultdict(float)
    # iterates through the kmer keys
    for kmer in variance:
        # deals with zero error division
        if variance[kmer] == 0.0:
            std[kmer] = std.get(kmer, 0.0)
        else:
            # calculates the standard deviantion and add 
            # the kmer and the standard deviantion values to the container            
            sd = math.sqrt(variance[kmer])
            std[kmer] = std.get(kmer, 0.0) + sd
    return std    
    
def z_scores(kmer_exp, kmer_counts, std):
    """
    Calculates the z scores to under/over represented kmers from a sequence.
    The score is calculaated as:
    
    Z(W) = (C(W) – E(C(W))) / sigma(W), where 
    C(w) - observed values
    E(C(w)) - represents the expected value from a kmer
    sigma - represents the standard deviation
    
    Inputs:
        kmer_exp - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    Outputs:
        z_scores - dictionary-like object mapping kmer to their z_scores.
    """
    # initialize the container
    z_scores = defaultdict(float)
    # iterates through the kmer keys
    for kmer in variance:
        # gets the kmer variance value
        var = variance[kmer]
        # gets the kmer std value
        sd = std[kmer]
        # deals with zero error division
        if sd == 0.0:
            z_scores[kmer] = z_scores.get(kmer, 0.0)
        else:
            # calculates the z score and add 
            # the kmer and the z score values to the container
            z = (kmer_counts[kmer] - kmer_exp[kmer]) / sd
            z_scores[kmer] = z
    return z_scores    
    
def get_p_values(z_scores_kmers):
    """
    Calculates the p value for all kmers.
    The calculation is done as:
    over represented: P(z > t) = erfc(t/sqrt(2))/2
    under represented: P(z > t) = erfc(-t/sqrt(2))/2
    t: thresholder
    
    Inputs:
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        
    Outputs:
        p_vals - dictionary-like object mapping kmer to their p values.
    """
    # initialize the container
    p_vals = defaultdict(float)
    # iterates through the kmer keys
    for kmer in z_scores_kmers:
        # calculates the p values to under represented
        # kmers (negative z scores)
        # add the kmer and p values to the container
        if z_scores_kmers[kmer] < 0.0:
            under = math.erfc(-z_scores_kmers[kmer] / math.sqrt(2)) / 2
            p_vals[kmer] = p_vals.get(kmer, 0.0) + under
        else:
            # add the kmer and p values to the container to over represented
            # and all other kmers
            others = math.erfc(z_scores_kmers[kmer] / math.sqrt(2)) / 2
    return p_vals    
    
def get_e_values(kmer_list, p_vals):
    """    
    Calculates the variance from a list of strings of length k.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        p_vals - a dictionary-like object mapping kmer to their calculated
                       p values.
    
    Outputs:
    
        e_values - a dictionary-like object mapping kmer to their calculated
                   e values.  
    """
    # number of tested hypoteses
    hyp_num = len(kmer_list)
    # initialize the container
    e_values = defaultdict(float)
    # iterates through the kmer list
    for kmer in kmer_list:
        # gets the p values from the input container
        p = hyp_num * p_vals[kmer]
        # calculates the e values and add the kmer 
        # and the e values to the container
        e_values[kmer] = e_values.get(kmer, 0.0) + p
    return e_values    
    
def gets_selected_kmers(kmer_list, 
                        kmer_counts, 
                        expected_kmers, 
                        z_scores_kmers, 
                        kmer_p_vals, 
                        kmer_e_vals, 
                        eval_cutoff=0.05):
    """
    Compile all the results from kmer analyses in a final csv for further analysis.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax. 
        kmer_exp - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.                   
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        p_vals - a dictionary-like object mapping kmer to their calculated
                       p values.
        e_values - a dictionary-like object mapping kmer to their calculated
                   e values.      
    Outputs:
    
        csv - a comma separated values with kmers and all calculated data.    
    """
    data = []
    to_check = []
    for kmer in kmer_list:
        keval = kmer_e_vals[kmer]
        if keval <= eval_cutoff:
            data.append((kmer, 
                         kmer_counts[kmer], 
                         expected_kmers[kmer], 
                         z_scores_kmers[kmer],
                         kmer_p_vals[kmer],
                         kmer_e_vals[kmer]))
        else:
            to_check.append((kmer, 
                             kmer_counts[kmer],
                             expected_kmers[kmer], 
                             z_scores_kmers[kmer],
                             kmer_p_vals[kmer],
                             kmer_e_vals[kmer]))
    df_filtered = pd.DataFrame(data, columns=['kmer', 
                                              'count',
                                              'expected',
                                              'z_score',
                                              'e_value',
                                              'p_value'
                                              ]).sort_values(by='z_score').reset_index(drop=True)
    df_tocheck = pd.DataFrame(to_check, columns=['kmer', 
                                                 'count',
                                                 'expected',
                                                 'z_score',
                                                 'e_value',
                                                 'p_value'
                                                 ]).sort_values(by='z_score').reset_index(drop=True)
    return df_filtered, df_tocheck    
    
def kmer_count(sequence, alphabet=['A', 'C', 'G', 'T'], k=1):
    """Return a dictionary with the k-mers as keys and it counts
    as values."""
    alphabet = set(alphabet)
    seq = sequence.upper()
    seq_len = len(seq)                                            
    kmers = [seq[i:i + k] for i in range(0, seq_len - k + 1)]     
    filterd_kmers = [kmer for kmer in kmers if                    
                     all(base in alphabet for base in kmer)] 
    return Counter(filterd_kmers)    
    
def get_count_frequencies(counts):
    """Receive a dictionary as key:counts."""
    total = sum(counts.values())
    return {key: (cnt/total) for key, cnt in counts.items()}    
    
def domainkey(m, symbol):
    # string module partition
    return m.partition(symbol)    
    
def variance(len_seq, prob_kmer):
    return len_seq * prob_kmer * (1 - prob_kmer)    
    
def expected_kmer_higher_markov(kmer_list, observed_count):
    """
    Calculates the expected value of a kmer based in a higher order Markov
    model. The model expand the middle for the two extremities of the kmer.
    
    Inputs:
    
        kmer_list - a list of kmers (string of length k).
        observed_count - a dictionary-likeobject with the observed/counts of
                         kmer (length between kmin and kmax) counted in a 
                         bigger sequence/genome.
    
    Outputs:
    
        expected_counts - a dictionary-likeobject with the expected counts of
                         kmer (length between kmin and kmax) according to a 
                         higher order Markov Model.        
    """
    # initialize the countainer
    expected_counts = defaultdict(int)
    # iterates through the content of the list
    for kmer in kmer_list:
        # gets the substrings
        pre, mid, suf = kmer[1:], kmer[1:-1], kmer[:-1]
        # look for the substrings counts and
        # get the expected values and add it to the countainer
        # it was add a small pseudo count to avoid erros by 
        # zero division
        expected = int((observed_count[suf] * observed_count[pre]) / (observed_count[mid] + 0.000001))
        expected_counts[kmer] = expected_counts.get(kmer, 0) + expected
    return expected_counts    
    
def variance_higher_order_markov(kmer_list, len_sequence, expected_count, k):
    """
    Calculates the variance of the HOMM as Var(C(W)) = N* P(W) * (1-P(W)) = E(C(W)) * (1 - E(C(W))/N),
    because the model is the sum of N almost independent observations, 
    each with probability P(W), it can be well modeled as a binomial distribution.
    
    Inputs:
    
        kmer_list - a list of kmers (string of length k).
        observed_count - a dictionary-likeobject with the observed/counts of
                         kmer (length between kmin and kmax) counted in a 
                         bigger sequence/genome.
    
    Outputs:
    
        expected_counts - a dictionary-likeobject with the expected counts of
                         kmer (length between kmin and kmax) according to a 
                         higher order Markov Model.        
    """
    std = defaultdict(float)
    for kmer in kmer_list:
        exp = expected_count[kmer]
        if exp == 0:
            sigma = 0.0
            std[kmer] = std.get(kmer, 0.0)
        else:
            sigma = math.sqrt(exp * (1 - exp/(N-k+1)))
            std[kmer] = std.get(kmer, 0.0) + sigma
    return std    
    
def get_z_scores(kmer_list, observed_count, expected_count, std):
    z_scores = defaultdict(float)
    for kmer in kmer_list:
        obs = observed_count[kmer]
        exp = expected_count[kmer]
        sigma = std[kmer]
        z_scores[kmer] = z_scores.get(kmer, 0.0) + (obs - exp)/sigma
    return z_scores    
    
def get_under_over_represented_kmers(kmer_list, len_sequence, z_scores, cutoff_e=0.05):
    evalues = defaultdict(float)
    random = defaultdict(float)
    for kmer in kmer_list:
        z = z_scores[kmer]
        under = len_sequence * math.erfc(-z/ math.sqrt(2)) / 2
        over = len_sequence * math.erfc(z/ math.sqrt(2)) / 2
        if under <= cutoff_e:
            evalues[kmer] = evalues.get(kmer, 0.0) + under
        elif over >= cutoff_e:
            evalues[kmer] = evalues.get(kmer, 0.0) + over
        else:
            evalues[kmer] = evalues.get(kmer, 0.0) + over
            evalues[kmer] = evalues.get(kmer, 0.0) + under
    return evalues, random    
    
def fasta_reader(filename):
    """
    Read a multi or single fasta file.
    
    Inputs:
    
        filename - string that represents a name of the file or a path to
                   the file.
    
    Outputs:
    
        A generator object containing a Seq and ID biopython objects.
    """
    from Bio.SeqIO.FastaIO import FastaIterator
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            for record in FastaIterator(handle):
                yield str(record.id), str(record.seq)
    else:
        with open(filename) as handle:
            for record in FastaIterator(handle):
                yield str(record.id), str(record.seq)    

def count_records_in_fasta_file(filename):
    if filename.endswith('gz'):
        with gzip.open(filename, 'rt') as fh:
            # alternate between a group of header lines
            # and a group of non-header lines
            # he number of groups is twice the number of sequences
            # len(list(itertools.groupby(f, key=isheader)))//2
            # but this is much better
            return sum(g for g,_ in itertools.groupby(fh, key=is_header))
    else:
        with open(filename, 'r') as fh:
            return sum(g for g,_ in itertools.groupby(fh, key=is_header))

def  kmers(kmin, kmax, seq): 
    """Get a 1st-order Markov model from a sequence of nucleotides."""
    for k in range(kmin, kmax + 1):
        for tup in tz.sliding_window(k, seq):
            yield ''.join(tup)

def get_sequence(path_to_files):
    """Stream a genome, letter by letter, from a list of FASTA filenames."""
    return tz.pipe(path_to_files,
                  cur.map(fasta_reader),
                  tz.concat,
                  cur.filter(is_sequence),
                  # concatenate characters from all lines
                   tz.concat,
                   # discard newlines and 'N'
                   cur.filter(is_nucleotide))
                   
def  genome_gz(file_pattern): 
    """Stream a genome, letter by letter, from a list of FASTA filenames.""" 
    return  tz.pipe(file_pattern,  
                    glob,  
                    sorted,   
                    # Filenames 
                    cur.map(gzopen(mode='rt')),   # lines 
                    # concatenate lines from all files: 
                    tz.concat, # drop header from each sequence 
                    cur.filter(is_sequence), 
                    # concatenate characters from all lines 
                    tz.concat, 
                    # discard newlines and 'N' 
                    cur.filter(is_nucleotide))                   
                   
# dicitonary of letters
LDICT = dict(zip('ACGT', range(4)))
LDICT                   
# dicitonary of of pairs of letters
PDICT = {(a, b): (LDICT[a], LDICT[b]) for a, b in it.product(LDICT, LDICT)}
PDICT                   
                   
def is_sequence(line):
    return not line.startswith('>')

def is_nucleotide(letter):
    # doesn't care abou Ns
    return letter in LDICT1                   
                   
@tz.curry
def increment_model(model, index):
    """  whenever  we  get  a  new  nucleotide  pair,  say,  (A, T), 
    we want to increment our Markov  model  (our  NumPy  matrix)  at  
    the  corresponding  position. """
    model[index] += 1                   
                   
def genome(file_pattern):
    """Stream a genome, letter by letter, from a list of FASTA filenames."""
    return tz.pipe(file_pattern, glob, sorted, # Filenames
                   c.map(open), # lines
                   # concatenate lines from all files:
                   tz.concat,
                   # drop header from each sequence
                   c.filter(is_sequence),
                   # concatenate characters from all lines
                   tz.concat,
                   # discard newlines and 'N'
                   c.filter(is_nucleotide))                   
                   
                   
                   

























































    
    
    
    
    
    
    
    
    
    
    
    
    
    

