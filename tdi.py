def get_observed_expected_ratio(kmer_count, expected, kmer_list):
    oe_ratio = defaultdict(float, [(k, 0.0) for k in kmer_list])
    for km in kmer_list:
        ob, exp = kmer_count[km], expected[km]
        if exp == 0:
            oe_ratio[km] = 0.0
        else:
            oe_ratio[km] = ob/exp
    return oe_ratio


def kmer_expected_zero_markov(base_freqs,  kmer_list, seq_len):
    """
    Calculates the 0-order markov given the background nucleotide distribution.
    
    Inputs:
        base_freqs - a dictionary-like obeject mapping the nucleotides to its calculated
                     frequence in a given sequence.
        kmer_list - a list of substrings of length k
        seq_len - a integer representing the length of the original sequence where the
                  calculated bases frequencies came from.
    
    Outputs:
        expected - a dictionary-like obeject mapping the expected kmer values
                   calculated as, (ex = math.pow(A, a)*math.pow(C, c)*math.pow(G, g)\
                   *math.pow(T, t)*(seq_len-k+1), where A,C,G,T are the frequency of each 
                   nucleotide in the genome, a,c,g,t are the number of each nucleotide in the k-mer,
                   and seq_len is the length of the kmer sequence.
    """
    # get kmer length
    k = len(kmer_list[0])
    # create the countainer and it keys from a kmer-list
    expected = defaultdict(float, [(km, 0.0) for km in kmer_list])
    # initiate the variables to receive the base frequencies
    A, C, G, T = base_freqs['A'], base_freqs['C'], base_freqs['G'], base_freqs['T']
    # iterate through each kmer in the list
    for kmer in kmer_list:
        # count the number of bases tha compound the kmer sequence
        a, c, g, t = kmer.count('A'), kmer.count('C'), kmer.count('G'), kmer.count('T')
        # calculates the expected values for taht kmer in the sequence accord to 
        # its base composition and add the value to the respectives keys in the dictionary
        ex = math.pow(A, a) * math.pow(C, c) * math.pow(G, g) * math.pow(T, t) * (seq_len - k + 1)
        expected[kmer] = ex
    # return the results
    return expected


def tud_normalized(kmer_count, expected_zero):
    """
    The normalized value for a word W is calculated by dividing the
    observed counts by the expected counts. This is the usage deviation
    vector for a genome. It is calculated as TUD(kmer) = Observed(kmer)/Expected(kmer).
    
    Input:
        kmer_count - a dictionary mapping the kmers found in a sequence to its counts (Observed).
        expected_zero - a dictionary mapping the kmers to its calculated expected values (Expected).
    
    Ouputs:
        tud_n - a dictionary mapping the kmers to its tetranucleotide usage deviation (TUD).
    """
    # create a list of kmer from the count dictionary
    kmer_list = list(kmer_count.keys())
    # create the keys and the countainer 
    tud_n = defaultdict(float,[(km, 0.0) for km in kmer_list])
    # for each kmer in the list
    for km in kmer_list:
        # get the expected and observed kmer values
        exp = expected_zero[km]
        obs = kmer_count[km]
        # to avoid zero division errors
        if exp == 0.0:
            tud_n[km] = 0.0
        else:
            # calculate the ratio of observed/expected values
            # add the values to it respective key
            tud_n[km] = obs/exp
    # return the result
    return tud_n


def TDI(sequence, oe_ratio, window, step, k):
    """
    Computes the k-mer difference index for a given genome. 
    Takes in a sequence, size of window to compute TDI in, astep size to move 
    along the genome by, k, a list of kmers of length k.
    
    Inputs:
        sequence - a string representing a DNA sequence
        oe_ratio - a dictionary-like object mapping the kmers to the
                   ratio observed/expected counts.
        
    """
	seq_len, start, end = len(sequence), 0, window
	w_i= []
	diff_windows = []
	while end < seq_len:
		if end > seq_len:
			end = seq_len
		seq_chunk = sequence[start:end]
		#record list of windows
		w_i.append([start,end])
		
		# calculate observed and expected in the window 
		# returns a dictionary 
		obs = kmerCount(seq_chunk, k)
		exp = zeroOrderExpected(seq_chunk, k)
		#compute difference sum for all kmers
		diff = 0
		for kmer in obs.keys():
			diff += abs((obs[kmer]/exp[kmer]) - oe_ratio[kmer])
		diff_windows.append(diff)

		#slide along the window by step
		start += step
		end += step

	xaxis = [x for x,y in w_i]
	return xaxis, diff_windows


def get_z_scores(diff_windows):
    # compute Z-score   Z=(x-mu)/sigma
	zscs = []
	mu = np.mean(diff_windows)
	sigma = np.std(diff_windows)
	for i in range(len(diff_windows)):
		zscs.append((diff_windows[i] - mu)/sigma)
		
		
def get_TDI(sequence, window, step, k, subset=False):
	if subset:
		sequence = sequence[subset[0]:subset[1]]    
	else:
	    tdi = TDI(sequence, window, step, k)

	return [tdi[0],tdi[1]]		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
