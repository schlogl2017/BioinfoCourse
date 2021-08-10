#!usr/bin/env python
# -*- coding: utf-8 -*-
#   python -di Data/Kmer_splitted -sd Chromosomes -ssd kmers \
#   -e csv.gz -do Results/kmers_merged

import os
import sys
import time
import glob
import argparse
from collections import defaultdict
from functools import reduce
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder, get_files, get_species_name
from system_utils import dask_csv_merge, concatenation_dfs

# def get_spc_df_merged(spc_name, filenames_dict):
#     dfs = []
#     for filename in filenames_dict[spc_name]:
#         dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
#     df_merged = pd.DataFrame(reduce(lambda left, right: pd.merge(left,
#                                                                  right,
#                                                                  on='kmer'),
#                                     dfs).set_index('kmer').sort_index().mean(axis=1).astype(int)).reset_index()
#     return df_merged.rename(columns={'kmer': 'kmers', 0: 'counts'}), len(dfs)


def parse_arguments():
    """Parse the command line arguments to the genome_palindromes_analysis script.
    Sets up the argparse command-line parser and calls it. These options can be
    accessed  using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to analyze palindromes and kmer in genomes.')
    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory of the files.')
    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='subdirectory of the files.')
    parser.add_argument('-ssd',
                        '--sub_sub_dir',
                        type=str,
                        dest='sub_sub_dir',
                        help='subdirectory of the files.')
    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-f',
                        '--text',
                        type=str,
                        dest='text_file',
                        help='Text file with genus names.')
    parser.add_argument('-e',
                        '--ext',
                        type=str,
                        dest='extension',
                        help='Final extension for the name of the files(i.e, .csv).')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    This script receive a path to the csv files, a directory name to save the
    merged csv files and a list of species names.
    The final result is a merged csv file with the mean of all kmer counts in the
    output directory.
    """
    print('Starting to process the script merge_csvs.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    # Data/Kmer_splitted (di)
    dir_in = opt.dir_in
    # Chromosomes (sd)
    sub_dir = opt.sub_dir
    # Results/Kmer_merged
    dir_out = opt.dir_out
    # kmers  (ssd)
    sub_sub_dir = opt.sub_sub_dir
    # text with names
    text_file = opt.text_file
    # Genus names
    names = get_species_name(text_file)
    # csv, csv.gz
    extension = opt.extension

    cnt = 0
    all_csvs = defaultdict(list)
    for name in names:
        print(f'Genus: {name}')
        all_csvs[name] = all_csvs.get(name, [])
        # Data/Kmer_splitted/Genus/Chromosomes/kmer
        full_name_dir = os.path.join(dir_in, name, sub_dir, sub_sub_dir)
        print(f'Getting csv files from species {name} from Directory: {full_name_dir}')
        for f in glob.glob(full_name_dir + f'/*.{extension}'):
            all_csvs[name].append(f)

    for name in names:
        # Results/Kmer_merged/genus/Chromosomes/kmer
        name_dir_out = os.path.join(dir_out, name, sub_dir, sub_sub_dir)
        # Genus_kmer_merged.csv
        csv_name = f'{name}_{sub_dir}_{sub_sub_dir}_merged.csv'
        # Results/Kmer_merged/genus/Chromosomes/kmer/Genus_Chromosome_kmer_merged.csv
        full_name = os.path.join(name_dir_out, csv_name)
        dfs, len_csvs = dask_csv_merge(name, all_csvs)
        # print(f'The number of {k}-kmers in these csv file is: {dfs["kmers"].size}')
        print(f'The number of csv files for species {name} is {len_csvs}\n')
        if not os.path.exists(name_dir_out):
            os.makedirs(name_dir_out)
        dfs.to_csv(f'{full_name}', index=False)
        print(f'The final csv merged file was saved in {full_name}\n')

        cnt += 1
    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {len(spc_names)}')
    print(f'Total files: {cnt}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!')


if __name__ == "__main__":
    sys.exit(main())
