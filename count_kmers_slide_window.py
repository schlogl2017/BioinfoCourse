#!/usr/bin/env python




def sequence_chunker(sequence, window_size, step):
    '''
    Return a iterator that yield start position, end position, and a property 
    calculated by a function, as tuple. The analysis is with overlapping windows.
    
    Inputs:
        sequence - a string representing a sequence (DNA/RNA/Protein)
        window_size - a integer representing the length of sequence to be analyzed
                   at each step in the full sequence.
        step - a integer representing the length of oerlapping sequence allowed.
    
    Outputs:
        A tuple:
            start - a number representing the start of the subsequence analyzed.
            end - a number representing the end of the subsequence analyzed.
            seq - string representing the subsequence/chunk
    '''
    seq = sequence.upper()
    
    for start in range(0, len(seq), step):
        end = start + window_size
        if end > len(seq):
            break
        yield start, end, seq[start:end]


def get_observed_expected_ratio(kmer_count, expected, kmer_list):
    oe_ratio = defaultdict(float, [(k, 0.0) for k in kmer_list])
    for km in kmer_list:
        ob, exp = kmer_count[km], expected[km]
        if exp == 0:
            oe_ratio[km] = 0.0
        else:
            oe_ratio[km] = ob/exp
    return oe_ratio


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
