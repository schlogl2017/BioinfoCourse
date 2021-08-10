#!usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
from time import time
import argparse
from termcolor import colored
from collections import defaultdict
from system_utils import make_me_a_folder, get_all_fasta
from fasta_parser import parse_fasta
from alphabet import iupac_dna
from sequence_utils import get_gc_content, count_all_bases, gc_content_sequence_window
from sequence_utils import base_content_slide_window
from system_utils import get_fasta_files


def parse_arguments():
    """Parse the command line arguments to the count_bases.
    Sets up the argparse command-line parser and calls it. These argsions can be accessed
    using args.argsion.
    The resulting results are csv files from each genome countained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(description="""A script to count all nucleotides from bacterial genomes/plasmids 
    in the fasta files.""", formatter_class=argparse.RawTextHelpFormatter)

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
    # alphabet
    alphabet = iupac_dna
    # get the list of all paths to the files in the input directory
    filenames = get_fasta_files(dir_in)
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
    for filename in filenames:
        print(colored(f"File: {filename}",
                      attrs=['bold']))
        # Data/Genomes_splitted/Genus
        # name of the taxon directory, ie. Acidisarcina
        # and get sub sub directory name
        genus = filename.split('/')[2]
        # read in the sequences and ids
        for seq_id, sequence in parse_fasta(filename):
            # get sequence length
            seq_len = len(sequence)
            print(f'Sequence length {seq_len}.')
            bases = count_all_bases(sequence)
            # Results/Genus/Bases
            path = os.path.join(dir_out, genus, sub_dir, sub_sub_dir)
            if not os.path.exists(path):
                os.makedirs(path)
            print(f'Saving the results in {path}\n')
            base_content_slide_window(sequence, path, seq_id, alphabet, 5000, 500, plot=True)
            with open(f'{path}/{seq_id}_bases.csv', 'w') as fout:
                fout.write('base,count\n')
                for base, cnt in bases.items():
                    fout.write(base + ',' + str(cnt) + '\n')
                    if not os.path.exists(path):
                        os.makedirs(path)
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
