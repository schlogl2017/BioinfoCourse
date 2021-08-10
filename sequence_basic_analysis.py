#!usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
import time
import itertools


start_time = time.process_time()


parser = argparse.ArgumentParser(description='Basic statistics on genomes.')
parser.add_argument('--path',
                    '-p',
                    required=True,
                    dest='path',
                    help='Path for the files')
parser.add_argument('--output',
                    '-o',
                    required=True,
                    dest='outputfile',
                    help='Name for output file')
parser.add_argument('--window_size',
                    '-ws',
                    type=int,
                    dest='window_size',
                    help='It represents the length of window size for genome analysis.')
parser.add_argument('--step',
                    '-s',
                    type=int,
                    dest='step',
                    help='Step for the window analysis.')
parser.add_argument('--GCcont',
                   '-gc',
                   action='store_true',
                   help='Option to return the GC content of a sequence.')
parser.add_argument('--ATcont',
                   '-at',
                   action='store_true',
                   help='Option to return the AT content of a sequence.')
parser.add_argument('--ATGCratio',
                   '-ratio',
                   action='store_true',
                   help='Option to return the AT/GC ratio of a sequence.')
args = parser.parse_args()


def is_header(line):
    return line[0] == '>'


def count_fasta_files(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        return sum(g for g, _ in itertools.groupby(f, key=is_header))


def str_punctuation_strip(word):
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for _ in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


def simple_fasta_reader(filename):
    print('Starting reading the fasta file' + '\n')
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as fh:
        fasta_iter = (it[1] for it in itertools.groupby(fh, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequence = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequence

def get_GC_content(sequence):
    """Returns the gc content of a sequence. It does a estimative
    of the GC content considering the putative N in tht total sequence."""
    seq = sequence.upper()
    len_seq = len(seq) - seq.count('N')
    gc = seq.count('G') + seq.count('C')
    return round((gc / len_seq) * 100, 4)


def get_AT_content(sequence):
    """Returns the gc content of a sequence."""
    seq = sequence.upper()
    len_seq = len(seq)
    at = seq.count('A') + seq.count('T')
    return round((at / len_seq) * 100, 4)


def get_ATGC_ratio(sequence):
    """It returns the GC/AT ratio from a sequence."""
    seq = sequence.upper()
    at = seq.count('A') + seq.count('T')
    gc = seq.count('G') + seq.count('C')
    return at / gc


def get_list_files(dirName):
    """Retruns a list of file and sub directories names in the directory"""
    list_files = os.listdir(dirName)
    all_files = list()
    # Iterate over all the entries
    for entry in list_files:
        # Create full path
        full_path = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files = all_files + get_list_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def check_filenames(files):
    zipped = []
    not_zip = []
    for file in files:
        if file.endswith('.gz'):
            zipped.append(file)
        elif file.endswith('.fa') or file.endswith('.faa'):
            not_zip.append(file)
    if not zipped:
        return not_zip
    return zipped


def sequenece_analysis_window_space(sequence, function, window_size, step):
    '''Return a iterator that yield start, end, and a property as tuple.
    The analysis is with overlapping windows.
    Ex.
    sequence = "attagcgcaatctaactacactactgccgcgcggcatatatttaaatata"
    for start, end, gc in sequenece_analysis_window_space(sequence, gc_content):
        print(start, end, gc)'''
    seq = sequence.upper()
    for start in range(0, len(seq), step):
        end = start + window_size
        if end > len(seq):
            break
        yield start, end, function(sequence[start:end])


# getting all the files
all_files = get_list_files(args.path)


print('Start reading the fasta files from directories: {}'.format(args.path))
files = check_filenames(all_files)


cnt = 0
for file in files:
    print('Fasta file: {}'.format(file))
    for name, seq in simple_fasta_reader(file):
        name_seq = name
        if args.GCcont:
            gc = get_GC_content(seq)
            analysis = sequenece_analysis_window_space(seq, get_GC_content,
                                                       args.window_size,
                                                       args.step)

        elif args.ATcont:
            at = get_AT_content(seq)
            analysis = sequenece_analysis_window_space(seq, get_AT_content,
                                                       args.window_size,
                                                       args.step)
        elif args.ATGCratio:
            atgc = get_ATGC_ratio(seq)
        cnt += 1


file_name = os.path.join(args.path, '.'.join((args.outputfile, 'tsv')))
with open(file_name, 'w') as fh:
    if args.GCcont:
        fh.write('start' + '\t' + 'end' + '\t' + 'data' + '\n')
        for star, end, data in analysis:
            fh.write(str(star) + '\t' +str(end) + '\t' + str(gc) + '\n')
    elif args.ATcont:
        fh.write('start' + '\t' + 'end' + '\t' + 'data' + '\n')
        for star, end, data in analysis:
            fh.write(str(star) + '\t' +str(end) + '\t' + str(at) + '\n')
    else:
        fh.write(f'The AT/GC ratio of the {name_seq} is: {atgc}')


print('Were read and manipulated: {} fasta files'.format(cnt))
print('Time of execution: {} seconds'.format(time.process_time() - start_time))
print('Script finished!')