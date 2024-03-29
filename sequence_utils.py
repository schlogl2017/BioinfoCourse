#!usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter
from itertools import product
import alphabet
from fasta_parser import fasta_item_counter, parse_fasta
from fasta_parser import parse_fasta
from count_kmers import count_kmers
from alphabet import iupac_dna

def get_chunks(sequence, window_size, step=1):
    """Returns a chunk of length of window_size and the end of the window size
    
    Inputs:
        sequence - a string representing a sequence.
        window_size - integer representing the length of the
                      chunk/window.
        step - a interger representing the length of the overlap window.
    
    Outputs:
        chunk - a string of length 'window'.
        end - a integer representing the end of the chunk.
    """
    # get the sequence length
    k = len(sequence)
    # get the index for each end and chunk
    for i in range(0, k - window_size + 1, step):
        # generate the end of the window
        end = i + window_size
        # get the slice of the sequence
        chunk = sequence[i:i + window_size]
        # assure the the chunk is the expected size
        assert len(chunk) == window_size
        yield chunk, end


def get_gc_content(sequence):
    """
    Finction to calculate the the gc content of a sequence.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.
    
    Outputs:
    
        gc - a float representing the of (g + c) content of a sequence.
    
    """
    # get the sequence length and 
    # make all the sequence characters upper case
    seq_len, sequence = len(sequence), sequence.upper()
    # count all gs and cs
    c = sequence.count('C')
    g = sequence.count('G')
    # returns the gc content from a sequence
    # sum up the |Gs and Cs counts and divide 
    # by the sequence length
    return round((c + g) / seq_len, 4)


def get_at_content(gc):
    """Returns the at content of a genome.
    
    Inputs:
    
        gc - a float representing the of (g + c) content of a sequence.
        
    Outputs:
    
        at content - returns the AT content of a sequence. It is
                     defined as 1 - the GC content.
    
    """
    return 1 - gc


def count_umbiguous_bases(sequence):
    """
    Function to count all umbigous bases in a sequence.
    Ambigous bases are bases that are not in the sequence
    alphabet, ie. 'ACGT' for DNA sequences.

    Inputs:

        sequence - a string representing a DNA sequence.

    Outputs:

        integer - represent the number of ambigous bases found and
                  counted in a sequence.
    """
    sequence = sequence.upper()
    amb = ['N', 'R', 'Y', 'W', 'S', 'K', 'M']
    return sum({base: sequence.count(base) for base in amb}.values())


def get_at_gc_ratio(at, gc):
    """
    Calculates the at/gc ratio of a genome.
    
    Input:
    
       gc - a float representing the of (g + c) content of a sequence.
       at - a float representing the of (a + t) content of a sequence. 
    
    Outputs:
    
       a float representing the ratio of the gc and at content of a sequence.
    
    """
    return at / gc


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


def gc_content_sequence_window(sequence, as_overlap=False, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. In
    overlapp windows of lenght k.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.    
        as_overlap - boolean that represents if overlap is needed.
        k - a integer reppresenting the lengths of overlappig bases.
            Default is 20.
    
    Outputs:
    
        gc_content - an array-like object with 
    
    
    """
    # make sequence upper case and getting the length of it
    sequence, seq_len = sequence.upper(), len(sequence)
    # the array-like object to collect the data
    gc_content = []
    # non overlap sequence length
    non_overlap = range(0, len(sequence) - k + 1, k)
    # overlap sequence length
    overlap = range(0, seq_len - k + 1)
    # overlap is needed
    if as_overlap:
        # iterates to the overlap region
        for i in overlap:
            # creates the substring to count the gc_content
            subseq = sequence[i:i + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    # if non overlap is choosed
    else:
        # iterates to the mon overlap region
        for j in non_overlap:
            # creates the substring to count the gc_content
            subseq = sequence[j:j + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    return gc_content


def codon_frequency(sequence, codon_table):
    """
    Function to calculate the frequency of the codons in a sequence.
    
    Inputs:
        
        sequence - a string representing a DNA sequence.
        codon_table - an array-like container with all possible codon, defined
                      as trinucleotides or triplets (64 possible triplets)
    
    Outputs:
    
        counter - a dictionary-like object with codon/triplets counts from a
                      given input sequence. 
        
    """
    # initialize the counter with the list of triplets from codon_table
    counter = Counter(dict([(c, 0) for c in codon_table]))
    # create a list/array of all possible codons found in the input sequence
    triplets = [sequence.upper()[i:i + 3] for
                i in range(0, len(sequence), 3)]
    # filters the triplets list from sequences that don't have length of 3
    # nucleotides
    triplets = filter(lambda x: len(x) == 3, triplets)
    # updates counter with the triplets counts and return it
    return counter + Counter(triplets)


def gc_var(sequence, as_overlap=False, k=20):
    """
    Calculates the gc content variance in a sequence according to a 
    window of length k.
    
    Inputs:
        sequence - a string representing a DNA sequence.
        k - integer representing the length of the search window.
            default is 20.
    
    Outputs:
    
        log of the gc variantion in the sequence in a window space of
        length k.
    
    """
    # calculates the percent of gc content
    gc = get_gc_content(sequence) * 100
    # get the gc content in the window space as an array
    gc_i = np.array(gc_content_sequence_window(sequence, as_overlap, k=k))
    # get the len of the gc content in the window space
    len_gc_i = np.shape(gc_i)[0]
    # check the difference of each point 
    dif = gc_i - gc
    return np.log((1 / len_gc_i) * sum(abs(dif)))


def base_stats(sequence, alphabet, as_count=False, as_dict=False):
    """Calculates de frequency or the number of bases in a sequence.
    
    Inputs:
    
        sequence - string representing the sequence
        alphabet - a alphabet (strings characters) that compound the string sequence
        as_count - boolean set as False
        as_dict - boolean set as False
    
    Output:
    
        counts - as default returns a numpy array as frequencies (floats) or
                 as a dictionary-like object
    
    Examples:
    
    > baseFreqs(seq, 'ACGT', asCounts = False, asDict = False)
    array([0.25, 0.25, 0.25, 0.25])

    as_count - True, returns a numpy array of counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = False)
    array([2, 2, 2, 2])

    as_dict - True and as_count as default (False) returns a dictionary as bases frequencies (float)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = False, asDict = True)
    {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    as_count True and as_dict True, returns a dictionary as base counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = True)
    {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    """
    # make the sequence upper case
    seq = sequence.upper()
    # count all bases in sequence and collect as an array
    counts = np.array([seq.count(i) for i in alphabet])
    # if is onle the counts
    if as_count:
        freqs = counts
    # other wise as frequencies
    else:
        freqs = counts / sum(counts * 1.0)
    # or as a dictionary like object
    if as_dict:
        return dict(zip(alphabet, freqs))
    else:
        return freqs


def get_strand_complement(sequence):
    """Returns the complement strand of the genome.
     
     Inputs:
        sequence - string representing the sequence   

    Outputs:
    
        sequence - string representing the complement of 
                   the string.    
    """
    # make the sequence upper case
    seq = sequence.upper()
    # table to change the complement characters
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """
    Returns the reverse complement strand of the genome.

    Inputs:

        sequence - string representing the sequence.

    Outputs:

        reversed_complement_sequence - string representing the reversed
                                       sequence complement.
    """
    return get_strand_complement(sequence)[::-1]


def get_sequence_skew(sequence):
    """
    Calculates the difference between the total number of
    occurrences of G and the total number of occurrences of C in
    the first i elements of the sequence. 

    Inputs:
        sequence - string representing the sequence     
    
    Outputs:
    
        skew - an array-like object that represents the GC skew of the sequence.
    
    Ex: 
    > get_sequence_skew('ACAACGTAGCAGTAGCAGTAGT')
    [0, 0, -1, -1, -1, -2, -1, -1, -1, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 2]
    
    """
    # make the sequence upper case
    sequence = sequence.upper()
    # start the array
    skew = [0]
    # iterates to the sequence elements and it indexes
    for idx, element in enumerate(sequence):
        # check if element[i] is a G
        # if so add 1
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        # if the element[i] is a C
        # add to the array -1
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            # if it is not G or C add 0
            skew.append(skew[idx])
    return skew


def get_minimum_skew(sequence):
    """
    Calculates a position in a sequence minimizing the skew.
    
    Inputs:
        sequence - string representing the sequence.    
    
    Outputs:
    
        min_skew - an array-like object that represents the pistion 
                   where the GC skew is the minimized in the sequence.    
    
    Example:
    get_minimum_skew('ACAACGTAGCAGTAGCAGTAGT')
    [5]
    """
    # start the array 
    min_skew = []
    # calculates the sequence gc skew
    skew = get_sequence_skew(sequence)
    # get the minimized skew values
    m_skew = min(skew)
    # iterates to the length of the sequence
    # to get the index positions
    for idx in range(len(sequence) + 1):
        # if the position i has the same value 
        # as the minimum appende to the array
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


def plot_base_frequency_genome(x_data, y_data, x_label, y_label):
    """
    Make a plot from the base frequency distribution in a DNA sequence.
    
    Inputs:
    
        x_data - sizes of the window where the frequencies were calculated
        y_data - base_freqs obtained in the windows of the sequence.
        x_label - the labels for the x axes
        y_label - the labels for the y axes
        
    Ouputs:
    
        plot of base composition in a sequence window.
        
    """
    # color for the bases
    base_markers = {"A": "b-",
                    "C": "r-",
                    "G": "g-",
                    "T": "y-",
                    "N": "k-"}
    # drawing the plot
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111)
    y_names = []
    for y in y_data:
        y_names.append(y)
        # adding colors to the lines representing the bases
        # the x and y data and the labels
        ax.plot(x_data, y_data[y], base_markers[y], label=y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    ax.legend(y_names)
    plt.grid(True)


def base_content_slide_window(sequence, path, name, alphabet, window, step, plot=False):
    """
    Calculates the base/nucleotide frequencies in a window of size window and
    step and make a plot of the base distribution along of the sequence length.
    
     Inputs:
     
        sequence - string representing the sequence.    
        name - if plot True is a string representing the name of the plot.
        alphabet - a alphabet (strings characters) that compound the string sequence
        window - integer representing the length of the search window.
        step - integer representing the size of sequence overlap in the window.
        plot - boolean value, default is False other wise True to draw the plot.
    
    Outputs:
    
        base_freqs - a dictionary-like object mapping the bases to its frequencies.
        sizes - array-like object tha represents the size of the windows
       
    """
    # sequence as a string of upper cases characters
    # bases as a set of upper cases characters
    sequence, bases = sequence.upper(), alphabet
    # initialize the dictionary container and the array
    base_freqs = defaultdict(list)
    sizes = []
    # iterates to the bases and start filling the dictionary
    # with the keys and a empty array
    for base in bases:
        base_freqs[base] = base_freqs.get(base, [])
    # iterates to the sequence windows
    for i in range(0, len(sequence) - window + 1, step):
        # gets the sequence of length of the desired window
        subseq = sequence[i:i + window]
        # check if the length of the window is correct
        assert (len(subseq) == window), 'The lenght of the subsequence must have the same size of the window'
        # start calculating the frequencies
        # and feeding the containers
        for base in bases:
            freq = subseq.count(base) / len(subseq) * 100
            base_freqs[base].append(round(freq, 4))
        sizes.append((i + window))
    # if it is to plot the data
    if plot:
        plot_base_frequency_genome(sizes, base_freqs, 'Genome (kb)', 'Frequencies')
        plt.title(f"Base Distribuition in {name} genome")
        plt.savefig(f"{path}/{name}_base_freq_slidewindow_plot.png")
    # return the data
    return base_freqs, sizes


def strand_stats(sequence, alphabet, start):
    """
    Calculates the DNA strand base statistics over a sequence.
    
      Inputs:
     
        sequence - string representing the sequence.    
        alphabet - a alphabet (strings characters) that compound the string sequence.
        start - integer representing the initial start position where to start
                calculate the bases statistics.
    
    Outputs:
    
        strand_stat - a dictionary-like object mapping the strands statistics.         
    
    """
    # assure the characters are upper case
    alphabet = alphabet.upper()
    # assure the characters are upper case
    # get the sequence length
    seq_len, seq = len(sequence), sequence.upper()
    # get the middle position of the sequence
    half_gen = (seq_len // 2)
    # get the final position
    ter = (start + half_gen)
    # initialyze the container
    strand_stat = defaultdict(tuple)
    # for circular genomes
    if ter > seq_len:
        ter = ter - (seq_len + 1)
    # iterates through the alphabet
    # count the bases
    for base in alphabet:
        base_total = seq.count(base)
        # check the strand
        if ter > start:
            for_strand = seq[start:ter].count(base)
            rev_strand = (base_total - for_strand)
        else:
            rev_strand = seq[ter:start].count(base)
            for_strand = (base_total - rev_strand)
        # calculates the differences between strand
        dif = (for_strand - rev_strand)
        # get the data in the container
        strand_stat[base] = (base_total, for_strand, rev_strand, dif)
    return strand_stat


def print_strand_stats(strand_statistics):
    """
    Prints the strand statistics.
    
      Inputs:
     
        strand_statistics - a dictionary-like object mapping the strands statistics.    
    
    Outputs:
    
        Prints out the strands statistics.    
    
    """
    print('   Total\tFor\tRev\tDif')
    for base, count in strand_statistics.items():
        print(f'{base}: {str(count[0])}\t{str(count[1])}\t{str(count[2])}\t{str(count[3])}')


def get_counts_from_kmer_list(filenames_lst, alphabet, kmin, kmax):
    """
    Prints the strand statistics.

      Inputs:

        filenames_lst - a array-like object that represents the paths to the
                        source of the fasta files.
        alphabet - a alphabet (strings characters) that compound the string sequence.
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        dic_list - a array-like object with kmers counts or each input fasta file.

    """
    # initialize the array container
    dic_list = []
    # iterates through the file paths
    for filename in filenames_lst:
        # get the sequences and ids
        for n, seq in parse_fasta(filename):
            # append the counts to the array
            dic_list.append(count_kmers(seq, alphabet, kmin, kmax))
    return dic_list


def get_mean_all_kmers_genomes_counts(filenames_lst, alphabet, kmin, kmax):
    """
    Count all kmers from a file source.

      Inputs:

        filenames_lst - a array-like object that represents the paths to the
                        source of the fasta files.
        alphabet - a alphabet (strings characters) that compound the string sequence.
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        Returns - a list of sorted tuples keys-values.

    """
    # initialyze the counter
    all_kmers = Counter()
    # get the number of files for each genus
    f_len = len(filenames_lst)
    # get the file path
    for filename in filenames_lst:
        # get the sequences
        for name, seq in parse_fasta(filename):
            # update the counter with kmer counts from the different
            # genomes in the directory source
            all_kmers.update(count_kmers(seq, alphabet, kmin, kmax))
    # get the average count from genomes in the directory source
    # for each genus
    kmer_all_counts = {k: (cnt // f_len) for (k, cnt) in all_kmers.items()}
    return sorted(kmer_all_counts.items(), key=lambda k: k[0])


def get_all_possible_kmers(alphabet, kmin, kmax):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet.

    Inputs:

        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmers - list of all possible combinations of k-mers of length k with length
                between kmin and kmax.

    """
    kmers = [''.join(letters) for n in range(kmin, kmax + 1)
             for letters in product(alphabet, repeat=n)]
    return kmers


def count_mers(sequence, alphabet, kmin, kmax):
    """Generate all DNA k-mers over the entirety of a sequence, with length
    between kmin and kmax.

    Inputs:

        sequence - string representing the sequence
        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmer_counts - a dictionary-like object with kmers counts from a
                      given input string with length between kmin and
                      kmax.

    """
    alphabet = set(alphabet)
    counts = defaultdict(int)
    for kmer in get_kmers_from_sequence(sequence, kmin, kmax):
        if set(kmer).issubset(alphabet):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def get_kmer_counts(kmer_list, kmer_counts):
    """
    Map all k-mers counts from a list of kmers of length k.

    Inputs:

        kmer_list - list or array-like object with kmers of length k.
        kmer_counts - a dictionary-like object with kmers counts from a
                      given input string with length between kmin and
                      kmax.

    Output:

        counts - a dictionary-like object with kmers counts length kma.

    counts will be utilyzed as input to calculate all kmer statistics.
    """
    counts = defaultdict(int)
    for kmer in kmer_list:
        counts[kmer] = counts.get(kmer, 0) + kmer_counts[kmer]
    return counts


def get_kmers_from_sequence(sequence, kmin, kmax):
    """
    Generate all DNA k-mers over the entirety of a sequence.

    Inputs:

        sequence - string where all kmers will be checked
        kmin: minimum DNA kmer length (int)
        kmax: maximum DNA kmer length (int)

    Output:

        yields all DNA kmers (str) of length kmin to kmax
    """
    limits = range(kmin, kmax + 1)
    seq_range = len(sequence) - kmax + 1
    for i in range(0, seq_range):
        for j in limits:
            yield sequence[i:i + j]


def get_pattern_positions(sequence, pattern):
    """
    Function to find a pattern in a string.
    
    Inputs:
        
        sequence - a string representing the sequence/genome.
        pattern - a substring to be found in the sequence.
        
    Outputs:
    
        tuple - a tuple with two integers representing the start and end 
                position of the pattern found in the sequence.
    OBS: re.finditer return an iterator yielding match objects over all 
    non-overlapping matches for the RE pattern in string. The string is 
    scanned left-to-right, and matches are returned in the order found.


    """
    return [pos.span()[0] for pos in re.finditer(r'(' + pattern + ')', sequence)]


def get_pattern_count(sequence, pattern):
    """
    Return the count the input pattern found in to a give string.
    :param sequence:
    :param pattern:
    :return:
    """
    return len(re.findall(r'(?=' + pattern + ')', sequence))


def get_kmer_count_slide_window(sequence, alphabet, window, step, kmin, kmax):
    slide_mer_count = defaultdict(Counter)
    for chunk, s, e in get_chunks(sequence, window, step):
        pos = '_'.join((str(s), str(e)))
        slide_mer_count[pos].update(count_kmers(chunk, alphabet, kmin, kmax))
    return slide_mer_count


def get_kmer_clumps(sequence, kmer_list, window, times):
    """
    Find clumps of repeated k-mers in string. A clump occurs when times or more 
    a k-mers appear within a window of size window. A list of (k-mer, position, count) 
    tuples is returned.
    
    Inputs:
        sequence - a string representing a sequence or a genome.
        alphabet - a set/list of all the allowed characters in the sequence.
        k - interger representing the kmers in the sequence. K-mers are 
            substrings of length k.
        window - a interger representing the length of the sequence to check by
                 the kmer clumps.
        times - a integer representig how many time each kmer must apear in the clump.
    
    Outputs:
         
         clumps - a dictionary-like object mapping the kmers to their positons and the end
                   of the clumps.
        
    """
    kmer_pos = defaultdict(list)
    k = len(kmer_list[0])
    clumps = defaultdict(list)
    for kmer in kmer_list:
        kmer_pos[kmer] = kmer_pos.get(kmer, []) + get_pattern_positions(sequence,
                                                                        kmer)
    for kmer, pos in kmer_pos.items():
        clumps[kmer] = clumps.get(kmer, [])
        for i in range(len(pos) - times):
            end = i + times - 1
            while (pos[end] - pos[i]) <= window - k:
                end += 1
                if end >= len(pos):
                    break
            if end - i >= times:
                clumps[kmer].append((pos[i], end - i))
    return clumps


def sequence_cleaner(sequence, alphabet):
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


def insert_wild_card(word, num_n=1):
    """
    Function to insert any wild card character in a word or string, generally
    the string must be a palindrome.

    Inputs:

        word - a string of characters that is best as a palindrome.
        num_n - a integer representing a the number of wild card characters
                to be inserted in the middle of the string word.

    Outputs:

        A tuple representing the word with n characters inserted in the middle and the
        inputed word without changes, if all definited condictions were achived, other
        wise return the word as it was inputed.
    """
    mid = len(word) // 2
    # to insert only one wild card character
    # with a predefinited condiction
    if num_n == 1 and is_palindrome(word) and len(word) % 2 != 0:
        return word[:mid] + 'N' + word[mid + 1:], word
    # the even words can receive two wild card chars
    # with a predefinited condiction
    elif num_n == 2 and is_palindrome(word) and len(word) % 2 == 0:
        return word[:mid - 1] + 'NN' + word[mid + 1:], word
    # only odd words can return a word with 3 chars wild cards
    # with a predefinited condiction
    elif num_n == 3 and word[:mid - 1] == get_reverse_complement(word[mid + 2:]) and len(word) % 2 != 0 and \
            len(word) >= 5:
        return word[:mid - 1] + 'NNN' + word[mid + 2:], word
    # if the condictions were not satisfied
    # return the word
    else:
        return word


def is_palindrome(string):
    """
     Function to check if a strings is palindromic or not.

     Inputs:

         string - a string of characters (a word, kmer, n-gram...)

     Outputs:

        boolean value - True if the string is palindromic other wise False.

     """
    k, mid = len(string), len(string) // 2
    # checking even palindromes
    if k % 2 == 0:
        return string[:mid] == get_reverse_complement(string[mid:])
    # checking odd palindromes
    else:
        return string[:mid] == get_reverse_complement(string[mid + 1:])


def get_counts(filename, alphabet, kmin, kmax):
    """
    Count all string countained in filename.

    Inputs:

        filename - a complete path to the file with th strings to count.
        alphabet - a alphabet (strings characters) that compound the string sequence
        min_k - minimum DNA kmer length (int)
        max_k - maximum DNA kmer length (int)

    Outputs:
        counter -  a dictionary-like object with kmers/string
                mapped to their counts

    """
    # get the list of kmers to count with length between kmin and kmax
    kmers_list = get_all_possible_kmers(alphabet, kmin, kmax)
    # initialyze the counter with all possible kmer with length
    # between kmin and kmax with zero counts
    counter = Counter(dict([(km, 0) for km in kmers_list]))
    # open and read in the kmers/string in the file
    with gzip.open(filename, 'rt') as fh:
        # iterates through the strings
        for line in fh:
            # make the adjustments int the strings
            kmer = line.replace('\n', '')
            # check if kmer/string is in the counter
            if kmer in counter:
                # if kmer is in add 1 other wise keep the zero count
                counter[kmer] += 1
    return counter


def get_counts_from_list(string_list, alphabet, kmin, kmax):
    """
    Count all string countained in filename.

    Inputs:

        filename - a complete path to the file with th strings to count.
        alphabet - a alphabet (strings characters) that compound the string sequence
        min_k - minimum DNA kmer length (int)
        max_k - maximum DNA kmer length (int)

    Outputs:
        counter -  a dictionary-like object with kmers/string
                mapped to their counts

    """
    # get the list of kmers to count with length between kmin and kmax
    kmers_list = get_all_possible_kmers(alphabet, kmin, kmax)
    # initialyze the counter with all possible kmer with length
    # between kmin and kmax with zero counts
    counter = Counter(dict([(km, 0) for km in kmers_list]))
    # open and read in the kmers/string in the file
    for string in string_list:
        # check if kmer/string is in the counter
        if string in counter:
            # if kmer is in add 1 other wise keep the zero count
            counter[string] += 1
    return counter


def get_mers(sequence, kmin, kmax):
    """returns the mers of length k obtained from a fasta file with genomes CDS or translated proteins."""
    for k in range(kmin, kmax + 1):
        return (''.join(mers) for mers in windowed(sequence, k))


def counts(sequence):
    """ 
    Find number of occurrences of each value in sequence.
    
    Inputs:
    
        sequence - a list or string like object.
    
    Outputs:
    
        count - a dictionary-like oblect mapping the items
                in the sequence to their counts (the frequency
                which they are fount in the sequence).
                
    >>> counts(['cat', 'cat', 'ox', 'pig', 'pig', 'cat'])
    {'cat': 3, 'ox': 1, 'pig': 2}
    
    >>> counts('cataa')
    {'c': 1, 'a': 3, 't': 1}
    """
    # initialize the countainer
    count = defaultdict(int)
    # iterates through sequence elements
    for item in sequence:
        # if element not in counts add 0
        # else add 1
        count[item] = count.get(item, 0) + 1
    return dict(count)


def count_subsequence_in_sliding_window(kmin, kmax, sequence):
    """ 
    Couns the number of overlapping subsequences of lenght between 
    kmin - kmax in a input sequence.
    
    Inputs:
    
        sequence - a string object representing a sequence.
        kmin - a integer representing the lower bound subsequence
               length.
        kmax - a integer representing the higher bound subsequence
               length.
               
    Outputs:
    
        Yields all subsequences of length between kmin and kmax found
        in a sequence.
    
    >>> list(count_subsequence_in_sliding_window(2, 3, 'ACGTACGTAACCCAGTGACC'))
    ['AC', 'CG', 'GT', 'TA', 'AC', 'CG', 'GT', 'TA', 'AA', 'AC', 'CC', 'CC', 'CA', 
    'AG', 'GT', 'TG', 'GA', 'AC', 'CC', 'ACG', 'CGT', 'GTA', 'TAC', 'ACG', 'CGT', '
    GTA', 'TAA', 'AAC', 'ACC', 'CCC', 'CCA', 'CAG', 'AGT', 'GTG', 'TGA', 'GAC', 'ACC']
    
    >>> ist(count_subsequence_in_sliding_window(2, 3, 'CCCAGTGACC'))
    ['CC', 'CC', 'CA', 'AG', 'GT', 'TG', 'GA', 'AC', 'CC', 'CCC', 'CCA', 'CAG', 'AGT', 
    'GTG', 'TGA', 'GAC', 'ACC']
    
    >>> list(count_subsequence_in_sliding_window(2, 3, ' ')) same as kmin, kmax = 0
    []
    
    >>> list(count_subsequence_in_sliding_window(0, 3, 'ACGT'))
    list(count_subsequence_in_sliding_window(0, 3, 'ACGT'))
    
    >>>list(count_subsequence_in_sliding_window(2, 2, ['1025']))
    []
    
    >>> list(count_subsequence_in_sliding_window(2, 2, '1025'))
    ['10', '02', '25']
    """
    if isinstance(sequence, str):
        for n in range(kmin, kmax + 1):
            for sub in zip(*(deque(itertools.islice(it, i), 0) or
                             it for i, it in enumerate(itertools.tee(sequence,
                                                                     n)))):
                yield ''.join(sub)


def get_gc_strand_difference(sequence, start):
    seq = sequence.upper()
    seq_len = len(seq)
    half = seq_len // 2
    ter = start + half
    g, c = 0, 0
    if ter > seq_len:
        ter = ter - seq_len + 1
    elif ter > star:
        g += 2 * seq[start:ter].count('G') - seq.count('G')
        c += 2 * seq[start:ter].count('C') - seq.count('C')
    else:
        g += seq.count('G') - 2 * seq[start:ter].count('G')
        c += seq.count('C') - 2 * seq[start:ter].count('C')
    return g - c


def plot_gc_skew(x, y):
    """Returns a plot of the genome gc skew
    input: GC_skew(sequence, n)"""
    plt.figure(num=None, figsize=(24, 7), dpi=100)
    yargmax = y.index(max(y))
    plt.axvline(oriCStart + oriOffset, color="r", linestyle='--')
    plt.axvline(x[yargmax], color="g", linestyle='--')
    plt.plot(x, y)


def get_neighbors(pattern, d):
    """
    Function to find a list of all offsets of patterns with hamming distance d of `pattern'.
    
    Inputs:
    
        pattern - substring represent a pattern
        d - integer representing the number of allowed charcaters that can be diferent
            from the pattern to they neighbors.
    
    Outputs:
    
        neighborhood - a sorted array-list-like object with all neighbors with d 
                       differences to the pattern.
    """
    # if no difference
    if d == 0:
        return [pattern]
    # if no pattern
    if len(pattern) == 1:
        return ['A', 'C', 'T', 'G']
    # initialize the container
    neighborhood = set()
    # checking for the suffix patterns
    neighbors = get_neighbors(pattern[1:], d)
    # iterates through the neighbors
    for kmer in neighbors:
        # check for the allowed distance
        if hamming_distance(pattern[1:], kmer) < d:
            # iterates through the charcater/bases
            for char in ['A', 'C', 'T', 'G']:
                # add the character to the suffix payyern
                neighborhood.add(char + kmer)
        else:
            # otherwise add the first character again
            neighborhood.add(pattern[0] + kmer)
    return sorted(list(neighborhood))


def count_sequence_mismatches(seq):
    """ When passed a string, representing a nucleotide sequence,
        treats it as a short inverted repeat, and returns the number of 
        mismatched compared to its reverse complement for half the length
        of the sequence.
    
    Inputs:
        sequence - a string representing a sequence.
        
    Outputs:
        mismatches - a integer representing the number of counted
                     differences between the two halfs of the sequence.    
    """
    trans_table = str.maketrans('ACGT', 'TGCA')
    half_len = len(seq) // 2
    second_half = seq[-half_len:].translate(trans_table)
    mismatches = 0
    for i in range(half_len):
        if seq[i] != second_half[-i - 1]:
            mismatches += 1
    return mismatches


def combination(n, k):
    """Combinations of N things taken k at a time.
    Exact long integer is computed.  Based on function from scipy.

    Notes:
      - If k > N, N < 0, or k < 0, then a 0 is returned.
    """
    if (k > n) or (n < 0) or (k < 0):
        return 0
    val = 1
    for j in range(min(k, N - k)):
        val = (val * (N - j)) // (j + 1)
    return val


def palindrome_search(sequence, min_len, max_len, alphabet, prob_cutoff=None):
    """
    Funtion to define the length of short inverted repeat (SIR) to looking for,
    with a min_len to max_len range of lengths.

    Inputs:
    
        sequence - a string representin a sequence of some length.
        min_len - a integer representing the minimum length of the 
                  palindrome.
        max_len - a integer representing the maximum length of the
                  palindrome.
        alphabet - a set/list/string representing all the characters
                   allowed in a sequence.Ex., 'ACGT' for DNA sequences.
        prob_cutoff - a float representing the probabilitie of that
                      paindrome occur given a background composition.
    
    Outputs:
    
        palindromes - a list/array with the length of the palindrome, 
        the start position, the probability of the palindrome, the
        number of mismatches and the palindrome sequence.
    """
    # get the sequence complement
    trans_table = str.maketrans('ACGT', 'TGCA')
    seq_complement = sequence.translate(trans_table)
    # gets the base composition
    nucs = base_stats(sequence, alphabet, False, True)
    # define maches bases
    matches = ['AT', 'TA', 'GC', 'CG']
    # probability of a match according tho the background
    p_match = 0
    # iterates tohrough the bases matches
    for b in matches:
        # calculate the probabilities
        p_match += nucs[b[0]] * nucs[b[1]]
    # checks if the results matches
    assert p_match == sum([nucs[b[0]] * nucs[b[1]] for b in matches])
    # initialize the container of possible probability using length and mismatches
    # as the indexes
    prob_dict = defaultdict(float)
    # iterates through the range of lengths
    for length in range(min_len, max_len):
        # iterates throught the half of the sequence
        for mismatches in range(0, (length // 2) + 1):
            # get the probabilities and number the mismatches
            p = probability(length, mismatches, p_match)
            prob_dict[(length, mismatches)] = prob_dict.get((length, mismatches), 0.0) + p
    # create an container for the results
    palindromes = []
    # iterates through the range of lengths
    for length in range(min_len, max_len):
        # defined mismatch threshold
        half_l = length // 2
        mismatches_cutoff = 0.5 * half_l
        half_list = range(half_l)
        # iterates throught to find the starts
        for start in range(0, (len(sequence) - length + 1)):
            # gets the putative palindromes
            seq = sequence[start:start + length]
            # gets the complement
            seq_comp = seq_complement[start:start + length]
            mismatches = 0
            # iterates throught the half lengths
            for i in half_list:
                # check for mismatches and increment the counts
                if seq[i] != seq_comp[-i - 1]:
                    mismatches += 1
            # check if the number of mismatches is allowed
            if mismatches <= mismatches_cutoff:
                # look up the probability,
                pr = prob_dict[(length, mismatches)]
                # if it passes the cutoff
                if pr <= prob_cutoff:
                    # add the results into the container
                    # count the number of the palindrome in the sequence
                    cnt_pal = get_pattern_count(sequence, seq)
                    palindromes += [[length, start, pr, mismatches, cnt_pal, seq]]
    return palindromes


def write_palindromes_to_file(path, csv_name, results):
    """ 
    Function to save the palindrome search result as a csv.
    
    Inputs:
        outfilename - the name of the final csv.
        results - a list of lists with the results of the palindrome
                  search. The data is the length, start position,
                  the probability of the palindrome, the mismatches number,
                  and the palindrome sequence.
    Outputs:
        A compressed or only a csv file.
    """
    if not os.path.exists(path):
        os.makedirs(path)
    df = pd.DataFrame(results, columns=['length',
                                        'start',
                                        'probability',
                                        'mismatches',
                                        'count',
                                        'sequence']).sort_values(by='probability').reset_index(drop=True)
    csv_name = f'{csv_name}_pal_search.csv'
    df.to_csv(f'{path}/{csv_name}.gz', index=False, compression='gzip')


def repeated_palindrome(palindromes_list):
    """
    Function to search through a data set (list of list) if one of the
    palindromes sequence has insert in it another palindrome.
    For example, if M is the longer palindrome sequence, if the p start > m start 
    and p end <  m end, p is inside m. The output file should be a new 
    data set containing only unique palindrome derived from the content of the input file
    
    Inputs:
        palindrome_list - a list containing the palindrome data as: 
        [[length, start, probability, mismatches, count, sequence]]
        
    Outputs:
        pal_list - a new palindrome list with unique palindromes as:
        [[length, start, probability, mismatches, count, sequence]]
    """
    # the list is ordered in the reversed form (long to short)
    ordered_palindrome = sorted(palindromes_list)
    longest_first = ordered_palindrome[::-1]
    # initialize a new list to receive unique plaindromes data
    pal_list = [longest_first[0]]
    # the longest palindrome cannot fit in any other sequence 
    # iterates over the longest_first original palindromes
    # get the start and end positions   
    for data in longest_first:
        start = data[1]
        end = start + data[0]
        # iterate through the pal_list and 
        # compare the start and end of the potential and palindromes 
        # to check if the potential palindrome is unique.
        unique_palindrome = None
        for dat in pal_list:
            start_unique = dat[1]
            end_unique = start_unique + dat[0]
            #  statement should test to check if the test palindrome fits
            # inside any of the identified 'real/unique' palindromes.
            if start >= start_unique and end <= end_unique:
                # if the palindrome tested fits inside
                unique_palindrome = False
                break
            else:
                # other wise it is unique
                unique_palindrome = True
        if unique_palindrome:
            # check if if it is not in the list
            if data not in pal_list:
                pal_list += [data]
    return pal_list


def check_the_palindromes_starts(palindromes_list):
    """ 
    This is a function to return start position minus start positions
    in an ordered data set of identified palindromes.
    """
    # to get all the start positions
    starts = [(s[1]) for s in palindromes_list]
    # all the data sorted
    palin = palindromes_list[1:]
    sorted(palin)
    # the sorted function orders the list for low to high
    longest_ordered = sorted(starts, reverse=True)
    # empty list to append the results to
    data = []
    # iteration through the start positions to return the results of the start
    # position minus the next start position in the list
    for i in range(len(longest_ordered) - 1):
        j = i + 1
        value = longest_ordered[i] - longest_ordered[j]
        data.append(value)
    return data


def get_cds_start_end_locations_genbank_file(filename):
    """
    Function to look for start and end of CDS feature in genbank record. 
    
    Inputs:
        filename - a path or locations of the genbank file.
    
    Outputs:
        genes - a dictionary like object that map CDS/genes to their
                locations in a genome.
    """
    # Loop over the features
    genes = defaultdict(list)
    cds = 0
    for seq_record in SeqIO.parse(filename, "genbank"):
        print(f'Dealing with GenBank record {seq_record.id}')
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS" and 'protein_id' in seq_feature.qualifiers:
                cds += 1
                prot_id = seq_feature.qualifiers['protein_id'][0]
                start, end = int(seq_feature.location.start), int(seq_feature.location.end)
                genes[prot_id] = genes.get(prot_id, []) + [start, end]
    print(f'There are {cds} CDS and {len(genes)} genes annoted for this genbank record')
    return genes


def final_kmer_counts(seq_dict, num_seqs, alphabet, min_k, max_k):
    """Function to count all the kmers found in a determinated sequence.
    
    Inputs:
    
        seq_dict - a dictionary-like object where keys are sequence ids and
                   values are the sequence itself.
                   
        num_seq - it a integer that represents the number of sequences that
                  are used to count the kmers
        
        alphabet - it is a string that represents the allowed characters
                   that compose the sequence. Ex: 'ACGT' for DNA sequences.
        
        min_k - minimum mer length (int)
        
        max_k - maximum mer length (int)
    
    Outputs:
        A tuple as:
         
             final_count - a dictionary-type with k-mers as keys and counts as values.
                           The final value is the average of counts from the total number
                           of sequences used to produce the counts.
             total_len - sum length is the average length of the sequences lengths 
    """
    counted = Counter()
    len_seqs = 0
    for name, sequence in seq_dict.items():
        seq = seq_cleaner(sequence, alphabet)
        len_seqs += len(seq)
        counted.update(count_kmers_cython(seq, min_k, max_k))
    final_count = {k: (v // num_seqs) for k, v in counted.items()}
    # total_len = (len_seqs // num_seqs)
    return final_count, len_seqs


def count_n_grams_fasta(fasta_dict, name, kmin, kmax):
    """
    Function to count all n-grams/k-mers (substrings of lenght n or k) in a
    big string/genome.

    Inputs:
        fasta_dict - a dictionary-like object that map a word/kmer to their value,
                    in this case a full path to the files to be analized.
        name - a string representing a word (key) that represent a key in a
               dictionary.
        kmin - a integer representing the lower bound of the kmer/n-gram length.
        kmax - a integer representing the maximum bound of the kmer/n-gram length.

    Outputs:
        final_counter - a dictionary-like mapping the kmers to their calculated count
                        in the input string, from a file.
    """
    # get the number of files in the names directory
    num_fastas = len(fasta_dict[name])
    # initialyze the counter
    counter = Counter()
    # iterates through the list of paths
    for filename in fasta_dict[name]:
        # reads the file and parse the content
        print(f'Reading and parsing the filename {filename}')
        for name, sequence in parse_fasta(filename):
            # counting the kmers
            cnt = count_kmers(sequence, kmin, kmax, counter=None)
            # add the count of the current file to the counter
            counter.update(cnt)
    # to get the mean of the kmer count for all the files
    final_counter = {k: (c // num_fastas) for k, c in counter.items()}
    return final_counter


def cleaning_ambiguous_bases(seq):
    """
    Function to clean up all umbigous bases in a sequence.
    Ambigous bases are bases that are not in the sequence
    alphabet, ie. 'ACGT' for DNA sequences.

    Inputs:

        sequence - a string representing a DNA sequence.

    Outputs:

        integer - a new clean up string representing a DNA sequence
                  without any ambiguous characteres.
    """
    # compile the regex with all ambiguous bases
    pat = re.compile(r'[NRYWXSKM]')
    # look for the ambiguous bases and replace by
    # nothing
    return re.sub(pat,  '',    seq)


def check_and_clean_sequence(sequence, alphabet):
    """
    Function to check and clean up all umbigous bases in a sequence.
    Ambigous bases are bases that are not in the sequence
    alphabet, ie. 'ACGT' for DNA sequences.

    Inputs:

        sequence - a string representing a DNA sequence.

    Outputs:

        cleaned_sequence - cleaned sequence or a string.
    """
    if set(sequence).issubset(alphabet):
        return sequence
    else:
        return cleaning_ambiguous_bases(sequence)


def cleaning_sequence_regex(sequence):
    """
    Function to check and clean up all umbigous bases in a sequence.
    Ambigous bases are bases that are not in the sequence
    alphabet, ie. 'ACGT' for DNA sequences.

    Inputs:

        sequence - a string representing a DNA sequence.

    Outputs:

        cleaned_sequence - cleaned sequence or a string.
    """
    amb = re.compile(r"[^ACGT]")
    return amb.sub("", sequence)

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


def do_bases_count_slide_window(sequence, window_size, step, alphabet, count=True):
    nuc_cnt_sw = defaultdict(dict)
    for st, end, seq in sequence_chunker(sequence, window_size, step):
        win = f'{st}-{end}'
        if count:
            nuc_cnt_sw[win] = base_stats(seq, alphabet, as_count=True, as_dict=True)
        else:
            nuc_cnt_sw[win] = base_stats(seq, alphabet, as_count=False, as_dict=True)
    return nuc_cnt_sw

def do_kmer_count_slide_window(sequence, window_size, step, k):
    km_cnt_sw = defaultdict(dict)
    for st, end, seq in sequence_chunker(sequence, window_size, step):
        win = f'{st}-{end}'
        km_cnt_sw[win] = count_kmers(seq, k)
    return km_cnt_sw
