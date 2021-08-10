#!/usr/bin/env python
# coding: utf-8
# usage: genome_kmer_statistics.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                                  [-ssd SUB_SUB_DIR] [-l CSV_FILENAME]
#                                  [-ka KMAX] [-ki KMIN] [-ev EVAL_CUTOFF]
# genome_kmer_statistics.py: error: the following arguments are required: -di/--dir_in
# python genome_kmer_statistics.py -di Results/Kmers_from_splitted -do Results/Kmer_statistics
# -sd Chromosomes -ssd kmers -l Results/Lengths/chr_lengths.csv -ki 5 -ka 5 -ev 0.01

import os
import sys
import argparse
from time import time
from termcolor import colored
from collections import defaultdict, Counter
from fasta_parser import parse_fasta, fasta_reader
from system_utils import get_len_csv, get_paths_to_csv_counts, \
    save_data_frame, make_me_a_folder
from sequence_utils import get_all_possible_kmers, get_kmer_count_from_csv, get_kmer_stats
from markov_models import get_expected_higher_markov, get_variance, get_standard_deviation
from markov_models import z_scores, get_p_values, get_e_values, gets_selected_kmers
from alphabet import iupac_dna


def save_data_frame(df, dir_out, name, sub_dir, sub_sub_dir, kmax, df_type):
    path = os.path.join(dir_out, name, sub_dir, sub_sub_dir)
    csv_name = f'{name}_{df_type}_{kmax}.csv'
    if not os.path.exists(path):
        os.makedirs(path)
    df.to_csv(f'{path}/{csv_name}', index=False)


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
                        help='Directory root. In my case the name is conjugated with a subdir, ex:Data/other.')

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
                        help='Name for a subdirectory, ex., kmers.')

    parser.add_argument('-l',
                        '--path_csv',
                        type=str,
                        dest='csv_filename',
                        help='Path to the csv file with sequence lengths.')

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

    parser.add_argument('-ev',
                        '--eval',
                        type=float,
                        default=0.01,
                        action='store',
                        dest='eval_cutoff',
                        help='Specify the e-value thresholder for select the k-mers.')

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
    # get the csv file with genome lengths
    csv_filename = args.csv_filename
    # cut off evalue
    eval_cutoff = args.eval_cutoff
    # check if the output directory existe other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)
    
    # Results/Lengths/chr_lengths.csv
    seq_len_dict = get_len_csv(csv_filename)
    # get the genus/species names
    names = seq_len_dict.keys()
    # get the csv files
    csv_files = get_paths_to_csv_counts(dir_in, sub_dir, sub_sub_dir, names)
    # get the kmer list
    kmer_list = get_all_possible_kmers(iupac_dna, kmin, kmax)
    # get all the stats and save it
    get_kmer_stats(seq_len_dict,
                   csv_files,
                   kmer_list,
                   dir_out,
                   sub_dir,
                   sub_sub_dir,
                   kmax,
                   eval_cutoff)
    # the final time
    end = time()
    # print some info
    print(colored(f'Total time for the script: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
