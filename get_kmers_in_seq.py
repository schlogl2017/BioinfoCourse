#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import time
from collections import Counter, defaultdict
from fasta_reader import parse_multi_fasta_file, str_punctuation_strip


parser = argparse.ArgumentParser(description='Counting and fetching the mers from genomes')
parser.add_argument('--path', '-p', type=str, required=True, dest='path',
                    help='Path for the fasta files')
parser.add_argument('--output', '-out', type=str, required=True, dest='output',
                    help='Filename for the csv output file')
parser.add_argument('--length', '-l', type=int, dest='length', required=True,
                    help='Integer representing the mers length')
parser.add_argument('--fasta', '-fa', action='store_true',
                   help='Return the list of translated cds from a genome.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--count', '-cnt', action='store_true',
                   help='If opted the script will count the mers.')
group.add_argument('--list', '-lt', action='store_true',
                   help='If opted the script will only returns a list of mer.')
args = parser.parse_args()


start_time = time.process_time()


def get_kmers(sequence, k):
    """Yields all the kmers of length k from the sequence input.
    sequence = a string representing a DNA, Protein or RNA
    k = integer representing the length of the kmer/substrings    
    """
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        yield kmer


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


# making a generator of the fasta files
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


dir_name = args.path
length = args.length


# getting all the files
all_files = get_list_files(dir_name)

print('Start reading the fasta files from directories: {}'.format(dir_name))
files = check_filenames(all_files)


mers_counter = Counter()
mers_lst = defaultdict()
cnt = 0
for file in files:
    print('Fasta file: {}'.format(file))
    for name, sequence in parse_multi_fasta_file(file):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        if args.count:
            mers_counter.update(list(mer for mer in get_kmers(sequence, args.length)))
        elif args.list:
            mers_lst[name] = ','.join(mer for mer in get_kmers(sequence, args.length))
    cnt += 1


output = args.output
print('Start writing the {} file.'.format(output))
with open(dir_name + '/' + output, 'w') as fout:
    if args.count:
        for (mers, count) in mers_counter.items():
            fout.write(mers + ',' + str(count) + '\n')
    elif args.list:
        for (name, mers) in mers_lst.items():
            fout.write(name + ',' + mers + '\n')
    elif args.fasta:
        for (name, mers) in mers_lst.items():
            fout.write(('>' + name + '\n' + mers + '\n'))


print('Were read and manipulated: {} fasta files'.format(cnt))
print('Time of execution: {} seconds'.format(time.process_time() - start_time))
print('Script finished!')