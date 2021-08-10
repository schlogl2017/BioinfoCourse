import re
import string

with open('dna.txt', 'rb') as f:
    data = f.read()

data = [x.split('\n', 1) for x in data.split('>')]
data = [(x[0], ''.join(x[1].split())) for x in data if len(x) == 2]

start, end = [re.compile(x) for x in 'ATG TAG|TGA|TAA'.split()]

revtrans = string.maketrans("ATGC","TACG")

def get_longest(starts, ends):
    ''' Simple brute-force for now.  Optimize later...
        Given a list of start locations and a list
        of end locations, return the longest valid
        string.  Returns tuple (length, start position)

        Assume starts and ends are sorted correctly
        from beginning to end of string.
    '''
    results = {}
    # Use smallest end that is bigger than each start
    ends.reverse()
    for start in starts:
        for end in ends:
            if end > start and (end - start) % 3 == 0:
                results[start] = end + 3
    results = [(end - start, start) for
               start, end in results.iteritems()]
    return max(results) if results else (0, 0)

def get_orfs(dna):
    ''' Returns length, header, forward/reverse indication,
        and longest match (corrected if reversed)
    '''
    header, seqf = dna
    seqr = seqf[::-1].translate(revtrans)
    def readgroup(seq, group):
        return list(x.start() for x in group.finditer(seq))
    f = get_longest(readgroup(seqf, start), readgroup(seqf, end))
    r = get_longest(readgroup(seqr, start), readgroup(seqr, end))
    (length, index), s, direction = max((f, seqf, 'forward'), (r, seqr, 'reverse'))
    return length, header, direction, s[index:index + length]

# Process entire file
all_orfs = [get_orfs(x) for x in data]

# Put in groups of 3
all_orfs = zip(all_orfs[::3], all_orfs[1::3], all_orfs[2::3])

# Process each group of 3
for x in all_orfs:
    x = max(x)   # Only pring longest in each group
    print(x)
    print('')
    
    
import io    # Only needed because input is in string form 
from Bio import Seq, SeqIO
import regex as re

startP = re.compile('ATG')

def get_orfs(nuc):
    orfs = []
    for m in startP.finditer(nuc, overlapped=True):
        pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
        orfs.append(nuc[m.start():m.start()+len(pro)*3+3])
    return orfs

for fasta in SeqIO.parse(io.StringIO(fasta_inputs), 'fasta'):
    header = fasta.description
    orfs = get_orfs(str(fasta.seq))
    print(header, orfs)   
    
"""Notes:

Normally you'd read a fasta collection from a file. Since here it's in string format, we used io.StringIO to make it easily compatible with SeqIO.parse from BioPython
The get_orfs function finds ATGs and returns the ORG originating from each one. If you're also interested in frames 4 through 6, you'll need the reverse_complement of the sequence.
If you're only interested in the longest ORF from each fasta sequence, you could have the get_orfs function return (max(orfs, key=len))
It's marginally more difficult if you're only interested in ORFs starting with ATG in a specific frame (e.g. frame 3). The simplest approach there might be to simply find all ORFs from the frame and then discard those not starting with ATG.
"""

import re
from string import maketrans

pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))')

def revcomp(dna_seq):
    return dna_seq[::-1].translate(maketrans("ATGC","TACG"))

def orfs(dna):
    return set(pattern.findall(dna) + pattern.findall(revcomp(dna)))

print orfs(Seq)


start = 'ATG'
stop = ('TAA', 'TGA', 'TAG')
stop_positions = []
start_positions = []

# Normal direction
for i in range(len(sequence)): # Create the index value
    if sequence[i:].startswith(start): # Check if there is a methionine in this position.
        start_positions.append(i)
    else:
        for stop_codon in stop:
            if sequence[i::].startswith(stop_codon): # Check if there is a stop in this position.
                stop_positions.append(i)
                break


for start_ in start_positions:
    for  stop_ in stop_positions:
        if stop_ <= start_:
            continue
        else:
            print("{}...{}".format(start_, stop_))
            break


import sys
from Bio import SeqIO
from Bio import Seq

record = SeqIO.read(open(sys.argv[1]), "fasta")
#Create three reading frames in forward direction, offset 0, 1, 2
readingframes = [Seq.translate(record.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]

results = []
for frame in readingframes:
    for peptide in frame.split('*'): #Split translation over stopcodons
        if len(peptide) > 30:
            results.append(peptide)

#Use PotentialORFs.txt as output, can be changed            
#Write length and translation to file
with open('PotentialORFs.txt', 'w') as output: 
    for peptide in results:
        output.write("{}\t{}\n".format(len(peptide), peptide))


def find_orf(sequence, gb):
        start_pos = sequence.find('GCCGCCACCATG')
        print "START " + str(start_pos)
        if start_pos >= 0:
                s_to_ATG = int(start_pos) + 9
                start = sequence[s_to_ATG:]
                for i in range(0, len(start), 3):
                        stops =["TAA", "TGA", "TAG"]
                        codon = start[i:i+3]
                        if codon in stops:
                                orf = start[:i+3]
                                break
                        else:
                                orf = start
                return orf, start_pos
        elif start_pos < 0:
                stop_pos = sequence.find(str(gb[-12:]))
                begin_to_stop = int(stop_pos) + 12
                return sequence[:begin_to_stop], start_pos

        else:
                print "Error: There is no open-reading frame for this sequence!"


# orffinder.py
# Author: Gungor Budak, www.gungorbudak.com
# Prints the longest ORF (in nucleotides) of a nucleotide sequence
# searching at six frames, where ORFs can start with any codon
# but ends with stop codons TAG, TGA, TAA. Six frames are the sequence
# itself and it's reverse complement, each of which has three frames
# by shifting one nucleotide right twice.
from string import maketrans
# Returns reverse complement of a sequence
def reverseComplement(seq):
	return seq[::-1].translate(maketrans("ATGC","TACG"))
# Marks stop codons searching the sequence in three nucleotides window
def markStopCodons(seq):
	edited_seq = ""
	i = 0
	while i < len(seq):
		# Codon
		codon = str(seq[i:i+3])
		# Replace stop codons with an asterisk
		if codon in ["TAG", "TGA", "TAA"]:
			edited_seq += "*"
		else:
			edited_seq += codon
		# Move to the next codon
		i += 3
	return edited_seq
# Return the longest sequence splitted from the nearest start codon to the beginning
def splitFromStartCodon(longest_orf):
	longest_orf_edited = ""
	i = 0
	is_start = False
	while i < len(longest_orf):
		codon = longest_orf[i:i+3]
		if codon == "ATG":
			is_start = True
		if is_start:
			longest_orf_edited += codon
		i += 3
	return longest_orf_edited
# Returns the longest ORF of a sequence which is also above the threshold
def findLongestOrf(seq, threshold):
	orfs = []
	for strand, nseq in [("+", seq), ("-", reverseComplement(seq))]:
		for frame in range(3):
			# Frame length to search codons in
			length = 3 * ((len(nseq)-frame) // 3)
			for orf in markStopCodons(nseq[frame:frame+length]).split("*"):
				if len(orf) >= threshold:
					orfs.append([orf, strand])
	max_orf = max([x[0]+x[1] for x in orfs], key=len)
	longest_orf, strand = max_orf[:-1], max_orf[-1]
	longest_orf_edited = splitFromStartCodon(longest_orf)
	return "Longest ORF: ", longest_orf_edited, "Strand: ", strand
# Main wrapper
def main():
	seq = "CATATTTAAAAATTAATTTACAAAGAAGACAGCGAGGCCGGGGTGGACGGCCGAGGTCAGGACGAGTCCGTGCAGGGGGACGGAGGGGGCGGGAAGAGGTGGCAGTAGCTTCTCTCAGTGACGGGCGAGGGCGGGCAGCTGTCCTCTGCAGACAGGTACTGGCTGTGGGGAGTGGGAGGGGGTGGGTAGGGGTCTGAGTCTGAATTCAAGTCCAGGTAGTATTTGCTGGTCTTCCAGCGACTGGTGCTGTAGTCACTGTCACACACGTCTGTGCTGCAAGGTGTCGTTGGGGGTGCCATGCCTCGGATGACGTAGGGCCTGCATGGTCTGGCGGTGGCAGGGATGTTTGAAGCGTAAAACATGTCCATGTTGTAGAGGGAGGGATCTGTGGCCGGGGATGGTGGTGGGTTCAGGATCGGTGGATACAGCGTGGCCTTTGTGCTGGATGAGCTGCTGGAGGAGGCCCCGGTGACGTGATTCCGGTCATAGAGGGGCACACTGCCCCGCCCCCCCACCAGGCTCATGGAATTCATCATGGACTTGCTGCACGGGATGCCTGGGAAGGGACCGTGCTGTGTGCCACCCGGGGCTATGAAGCTGAGAGGCACATGAGGAGTGCCACTGACATACTCATGGGGAAAGGGCCCGTTGGCTCCCGTGTAGCGCTGGCACATCACACGCTGGCAGACAAAGTAGACCCCGCCCATGACGAAGAGGGAGAGGATGATGCCGATGACAGGCCCGATGGCGCTGCTGTGGGCTGGCGTGTCATCAGAGGGCGGCTTGTTGATTTCACACATGAGCTCGTCAGATCCGTCGGCGCAGTCGGGGAAGGAGTCGCACTGCTGCTTGATGTGGACGCACTGGCCGCTGGCACACCGGAACTGATTGGGCAGACAGACAGCGTCGCAGTTAGCTTCATCGGAGCGGTCCTGGCAGTCGGCTTCGCCGTCGCAGCGCAGCCGCAGGTCCACGCACTGGCCCCGAGCGCAGGGGAACTGGGAGGCAGAGCACACTGGGCAGCCCTCCTCATCGCTCTGGTCTGCACACTCAGGGAAGCCGTCACAGCGCCAGGCTCCAGGGATGCAGTCAATCTCACCGGTGGCACAAGCAAACTGATCAGGGGAGCAGGTAGGAGGCTCGCCACAAGTCAGCAGGTTCTGCAAGAGCACCAGGTGGACTGGGCAGGAGCAACGCGGCGTTCCATCCCCCTTGGCGATGCAGATGTGGGAGCAGCCTCCATTGTCTCGGGCACACGGATGGGCTGAGAATTCCTCCAGGCTGACTTCCTCCACGGCATGGATGCCGGTCAGATGGGTGACACGGCCCTGAACCCTAGTCCGCTTGTCCCCGGTGGTCTTCTCCACCCGCTCGATCATCTGCTGCTGCCGGTCGATCCAGTAGAGGTGCCTGCCCAGCACCGTCAGGCCCACTGGCTGCACGATGTTGGCATCTTCCAGGGTCAGACGGTTGGCCCCCGACAGGTCGCAGCTTTCGATGCGCTTTAGGTCTGCGTCCACCCAGAACAGCTTGCCCAGAGCATTGTCCACCACAAGGGCCACAGGACGGATGAGGCCCGTGGTGAAGAGAACCTCACGCTCCGTGCCATCCAGGGAGGCCCGCTCGATCTTGGCCGCGTGGTCCTGCATATTGGTGAAGTACATGTACCCTCGCTCAGCATTGACAGCAATGGCCCTTGGCTTGTCACGGTCCCCTCGAAGAACCACTCCCATGGCATCCCCATCCAGCCGGTGGACATTGATGGTATTGGTGGCCTCACAGGTCCAGAACAGTGTCCGGCTGTAGATGTCAATGCTGAGATCGTGTGGCTGTCTGTCTGGGCTTAGGCCTTGGCTGGGAGAGGTCAGCACGGAGGGCTGGGTCCCGTCGTCCTTGGCCCTCTTGATGTTCTGGCGCCCGTCCACCCAGTAGATGAACTTGTCCAGCGGGTCATAGTTGATGGCTTTGACGTTCCTCAGCCCGTGCAGGGGCAGGACAAGGTCCGGGCTGAGCTGGTCATCGGGGATCATCCGGCTGATGGCAAATTTCTGGCTGAACAGCAAGAAGGTGGAGGGCGGGCTACAGTTGCGGCTGCTGGGGTCCAGCGTGTAGTGAGAAGCACAGCCACAGCGGTGGCCCCCGGGGATGGCGAGGCACAGCTGCCCGCACTGGCCATTATTGTGCACGCAGTCGTTGAGGCCATCCTGACGGGAGGAGTGGAACACCAGAAGGTCCATGACGAAGTCCAGGTGGCCCTGGATGAGGGTGCGGTTCCGCCCGCTGGTCTTGTCGGCCCGCTCGATACTGTGCAGGTTCCAGTCAGTCCAGTAGATGTAATCGCTGTACTGAGTCAGGCCAAAGGGGTAGGGCAGATCGTCGGCTATCACCATGCGCTCCTGACCCAGCATGTTAGAGGACTCAATCATGTTGGTGTCCAGGTCAGTCCAGTACAGTCGCTGGTCGGCGTAGTCGATGGTGAGGTCATTGGCCCGGCCCACCTTGTCTACCAGCGTCATGCAATTGGTCCCGTCCATGAAGGCCCGCACAATCCTTGGCTTGCCACCCCACTCGGTCCAGTAGATGTAGCCTTTAGTAGGATCCAGAGCCAGAGACCTGGGGTTGTCGAGGTCTCTCCACACGAGCACTTGGCGGAACTGCCCATCCAGCCGGGCCACCTCAATCCTGTTGGTCCCTGTGTCAGCCCAATAGAGGTTCTTGCCCATCCAGTCCACGGCCATTCCCTCAGGGTAGTCGAGACCAAACTCGATCACGTGCTCCACTGAGCTCCCATTCATGAAGGCTCGGCTGATAGTCTTGAGGCTGACATCGGTCCAGTAGATGTGGTTGTTGGACACATCAAAGTCCAGGGCAGAGGCCTCCTTGACACCTGTGAGTGGGATGGCCACATCGTTGTTGTTAGTCTCCAGGGAGATCCTGTGGATGGTGGCTCTGCTGGTGAACACCAGGAAGGCCTCGGGGACAATGCAGGTCTTCATGTCACTCAACAGCTCCAGGCCGATGGGGCAGCCACACTTGGTGGCATGAGGGGTGAAGAAACACAGATGGCTGCACCCTCCGTTCCCGTCTGCGCATGGGTTGGTTCCAACAACCTTGACCACATTCACTGCTTTGAGTCCCATCAGGTCGGGCAGCTGATCGATGATGACATCCCGGCTGGCCTTGACCTTGTGTACCCTCTCAATACTGCGCCGCTGCCAGTCGGTCCAGTAGATGAAGTCCCCCAGCAGTGTGAACCCGAAAATGTGTGGGAGCTTGTCTTCTAGCAGGGTCCGCCGCTTTGTCCCGTCTATGTTGATCACCTCGATTTTATCAGTCTTGGCATCCCCCCAGTAGAGCTTGCCTTCCTGCAGATCCAGGGCCAGCCCGTTGGGCCACCCAAGGGAGGTGTTCACAAGGACATGCCGATCTCGTCCATCTAGGTTGGCACATTCGATTTTGGGGTTCTCCCCCCAATCTGTCCAGTACATGAGGCCCATCACAGGGTGCAACACAATGGCCCGAGGTTCGTCCAGGTCCTCAGATACCAGGATCTTTCGGGAGGTACCGTTGAGGCGAGTCACCTCAATTCGGTCAGTGCCTGTATCTGTCCAGTAGAGGTTCCGGGCGACCCAGTCCACAGCAATACCATCGGGGTCATTGATCTCAGTGTTCACAAGTGTCTGTGCGCCTGACCCATCCAGGTATGCCCTGCGGATAGCCCGCACCTCATCGTCCGTCCAGTACACGTAGCCCTCCAGCGGGTCGTAGTCAATGGCAATGGCATGCCGGATGTCTCCCACCTGCAGCACTATGTCTGTGAAGTCAGGTGTGTCCAAAGAGATCCTCCTCAGGTCTGTCCTCCGGGCCAGCAGCAGCACTTCCTCAGCCCCTTCCTTGCACGTCTTGCCATTGTCCCGCAACTGGACACCAGTAGGGCAGGCGCAGGAGTAGAACGGCTCCCGTGGGGACAGCAGGCATAGGTGGGAGCAGCCGCCGTTGTCCTCCTCGCATCGCGTGTGGAAGGGAGGCTGGCGCTCCTGGCTCAGCACCTGAATGTCCATGGGTGAGTAGAGGGCACTGAGGATCTCCTTCCTCTTCTCCCCTGTCCACTTGTTGCAGGCGTGGATAGAGCGGGTCTGCCAGTCCGTCCAGTAGAGTGTGTCGCCAGAGAGCGTCAGGGCAAAGGGGTGAGTGAGGCTGCCCTCCACCACCTTCTGCCGGAAGGAGCCGTCCAGGTTGGCACGGTGGATGAAGCTGAGCTTGGCATCAGCCCAGTACAGCTTCTGCTCCTCCAGATCGATGGTCAGTCCGTTGGGCCAGTAAATGTCAGAGTCCACGATGATCTTCCGGGTACTGCCATCCATCCCCGCCCGCTCGATCCGGGGTGCCTCCCCCCAGTCAGTCCAGTACATGTACCCGTGTGCAGGATCCAGGGCAATGGCCCTTGGCTGGTCCAGGTCTTGCCAGAAGAGAACCTTACGGGACGTCCCATTGAGGTTGGCCACCTCGATGCGGTTGGTCTCTGAGTCTGTCCAGTACAGCTTCTTGCCAACCCAGTCACAGGCCAGGCCATCGGGTGACACGAGGCCGGAGATGACAATGCTCTGCGTGGCAGTTCCATTCTGGCTTAGGAAGGTCTGTTTGATGGCCTCCTCGCTCACGTCCGTCCAGTACACAGCGCCCTTGGAGAACTGGAAGTCCACGGCAGCCGCATCCTCCAGGCCACTGGCTACAATGGTGGCCTCCAGCTTCACCCCACCGGCATCCACTAGCCGCACATCCCGGCGGTTGGCAAACAACAGGAGCGGTGAGGCCGCGGCGGGGGCCAGGCTGCAGCACAGCAACAGCAGCAACAGCAGCGGCGGCGGCAGCGGCGGCAGAGGGGCTCGGGTCGGCACCGTTTCCATGTCGTCCGGCGGCCGGGAGCGCCGCGGGCTCACTCGGGCTCCATGGCGCGCGCGGCGGCGGCAGCGGCGACGGCGGTGGATCCTCGCGGCTCCCGGCGCCTCCTGCTCCCGCCTCGGG"
	print findLongestOrf(seq, 180)
# Run
main()









 
    
