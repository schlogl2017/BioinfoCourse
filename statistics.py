#!usr/bin/env python
# -*- coding: utf-8 -*-
# usage: statistics.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                      [-ssd SUB_SUB_DIR] [-e EXTENSION] [-ka KMAX] [-ki KMIN]
# statistics.py: error: the following arguments are required: -di/--dir_in
# python statistics.py -di Data/Genomes_splitted -do Results/Kmer_statistics -sd Chromosomes -ssd kmers -e gz
# -ki 6 -ka 6


import os
import sys
import gzip
from time import time
import pandas as pd
from termcolor import colored
import argparse
from collections import defaultdict, Counter
from system_utils import make_me_a_folder, get_all_fasta, save_data_frame_kmers
from fasta_parser import parse_fasta
from count_kmers import count_kmers
from sequence_utils import count_k_mers_fasta, split_seq, get_all_possible_kmers
from markov_models import get_expected_higher_markov, get_variance, get_standard_deviation, \
    z_scores, get_p_values, get_e_values, get_kmer_data
from alphabet import iupac_dna


def parse_arguments():
    """Parse the command line arguments to the the statistics.
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
                        '--extension',
                        type=str,
                        dest='extension',
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
    # extension type for fasta
    extension = args.extension
    # alphabet
    alphabet = iupac_dna
    # get the list of all paths to the files in the input directory
    # ex., Data/Genomes_splitted, Chromosomes, gz
    fasta_dict = get_all_fasta(dir_in, sub_dir, extension)
    # check if the output directory existe other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # initialize the file counter
    cnt_files = 0
    # input the file paths and print it to show where the script is doing
    for name in fasta_dict.keys():
        print(colored(f"Start working with genus {name}\n", attrs=['bold']))
        # initialize the kmers counts
        cnt, seq_len = count_k_mers_fasta(fasta_dict, name, alphabet, kmax-2, kmax, overlap=kmax, nprocs=4)
        # get the k-mer list for analyziz, k = 6
        kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)
        # calculating the expected number for all k-mers
        expected = get_expected_higher_markov(kmer_list, cnt)
        # get the expected count variance
        variance = get_variance(kmer_list, seq_len, expected)
        # get the standard deviation
        std = get_standard_deviation(variance)
        # getting the z-scores
        z_scrs = z_scores(expected, cnt, std)
        # get the p-values from k-mers
        pvals = get_p_values(z_scrs)
        # get the k-mers e-values
        evals = get_e_values(kmer_list, pvals)
        # saving the final results as a csv file
        kmers = get_kmer_data(kmer_list,
                              cnt,
                              expected,
                              z_scrs,
                              evals,
                              pvals)
        save_data_frame_kmers(dir_out, sub_dir, sub_sub_dir, name, kmax, kmers)
        print(f'Number of kmer (kmin-{kmax-2}/kmax-{kmax}) from {name}: {len(cnt)}\n')
        # k = kmax
        # k_mers = len(expected)
        # pos = seq_len - k + 1
        # all_mers = 4 ** k
        # mis = (4 ** k) - k_mers
        # rep = ((seq_len - 1) - k + 1) - k_mers
        # with open(f'Results/{name}_kmers{kmax}.txt', 'w') as fh:
        #     fh.write('k\tkmers\t4^k\tpositions\tmissing\trepeated\n')
        #     # (k, len(kmers), 4**k, (len(seq[0])-1)-k+1, 4**k-len(kmers), (len(seq[0])-1)-k+1-len(kmers))
        #     fh.write(f'{k}\t{k_mers}\t{all_mers}\t{pos}\t{mis}\t{rep}\n')
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
