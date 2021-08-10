#!usr/bin/env python
# -*- coding: utf-8 -*-

# usage: count_all_kmers_from_fasta.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                                      [-ssd SUB_SUB_DIR] [-e EXTENTION]
#                                      [-ka KMAX] [-ki KMIN]
# count_all_kmers_from_fasta.py: the following arguments are required: -di/--dir_in
# dir_in, subdir, extention
# get_all_fasta('Data/Genomes_splitted', 'Chromosomes', 'gz')


import os
import sys
import gzip
from time import time
import pandas as pd
from termcolor import colored
import argparse
from collections import defaultdict, Counter
from system_utils import make_me_a_folder, get_all_fasta
from sequence_utils import count_k_mers_fasta
from fasta_parser import parse_fasta
from count_kmers import count_kmers
from alphabet import iupac_dna


def count_n_grams_fasta(fasta_dict, name, alphabet, kmin, kmax):
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
    # alphabet as a set
    alphabet = set(alphabet)
    # get the number of files in the names directory
    num_fastas = len(fasta_dict[name])
    print(f'The number of fasta files for this genus is {num_fastas}.')
    # initialyze the counter
    counter = Counter()
    # iterates through the list of paths
    for filename in fasta_dict[name]:
        # reads the file and parse the content
        print(f'Reading and parsing the file {filename}')
        for name, sequence in parse_fasta(filename):
            print(f'Sequence length {len(sequence)}')
            # get the counting the kmers
            cnt = count_kmers(sequence, kmin, kmax, counter=None)
            # add the count of the current file to the counter
            counter.update(cnt)
    # to get the mean of the kmer count for all the files
    final_counter = {k: (c // num_fastas) for k, c in counter.items() if set(k).issubset(alphabet)}
    return final_counter


def parse_arguments():
    """Parse the command line arguments to the the count_all_kmers_from_fasta.
    Sets up the argparse command-line parser and calls it. These argsions can be accessed
    using args.argsion.
    The resulting results are csv files from each genome countained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(description="""A script to generate all possible kmers
    of length k from bacterial genomes/plasmids in the fasta files.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory root. In my case the name is conjugated with a subdir')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Name for a subdirectory, ex., Chromosomes.')

    parser.add_argument('-ssd',
                        '--sub_sub_dir',
                        type=str,
                        dest='sub_sub_dir',
                        help='Name for a subdirectory, ex., Chromosomes.')

    parser.add_argument('-e',
                        '--extention',
                        type=str,
                        dest='extention',
                        help='Name representing the type file. Ex., gz')

    parser.add_argument('-ka',
                        '--kmax',
                        type=int,
                        default=6,
                        action='store',
                        dest='kmax',
                        help='Specify the maximum kmer/palindrome length.')

    parser.add_argument('-ki',
                        '--kmin',
                        type=int,
                        default=6,
                        action='store',
                        dest='kmin',
                        help='Specify the minimum kmer/palindrome length.')

    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input diretory, ex. Data/Genomes_splitted
    dir_in = args.dir_in
    # name of the sub directory to save the final result
    # Chromosomes/Plasmids
    sub_dir = args.sub_dir
    # sub_sub dir name, ex., kmers/palindromes
    sub_sub_dir = args.sub_sub_dir
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # minimum kmer length
    kmin = args.kmin
    # miximum kmer length
    kmax = args.kmax
    # extention type for fasta
    extention = args.extention
    # alphabet
    alphabet = iupac_dna
    # get the list of all paths to the files in the input directory
    # ex., Data/Genomes_splitted, Chromosomes, gz
    fasta_dict = get_all_fasta(dir_in, sub_dir, extention)
    # check if the output directory existe other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # initialyze the file counter
    cnt_files = 0
    # input the file paths and print it to show where the script is doing
    for name in fasta_dict.keys():
        print(colored(f"Start working with genus {name}\n", attrs=['bold']))
        # initialyze the kmers counts
        cnt = count_n_grams_fasta(fasta_dict, name, alphabet, kmin, kmax)
        print(f'Number of kmer (kmin-{kmin}/kmax-{kmax}) from {name}: {len(cnt)}.\n')
        # get the name of the full path to save the final csv file
        # Results/genus/Chromosomes/kmers{k}/ids{k}.csv
        full_path = os.path.join(dir_out,
                                 name,
                                 sub_dir,
                                 sub_sub_dir)
        # checking if there are a path to save the data
        # if not create it
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        # name of the file to be saved
        csv_name = f'{name}_{sub_dir}{kmax}.csv'
        # iterate through the list and write the kmer to the file
        df = pd.DataFrame(cnt.items(), columns=['kmer', 'count'])
        df.to_csv(f'{full_path}/{csv_name}.gz', index=False, compression='gzip')
        # add to the count file
        cnt_files += 1
    # the final time
    end = time()
    # print some info
    print(colored(f"Total number of files: {cnt_files}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
