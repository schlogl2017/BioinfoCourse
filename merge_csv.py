#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import glob
import argparse
from collections import defaultdict
from functools import reduce
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder, get_files

def get_species_name(filenames):
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def parse_arguments():
    """Parse the command line arguments to the genome_palindromes_analysis script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to merge csv files.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Directory of the files.')
    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='bsudirectory of the files.')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-e',
                        '--ext',
                        type=str,
                        dest='extension',
                        help='Final extention for the name of the files(i.e, .csv).')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        dest='filenames',
                        help='species list  as txt.')
    parser.add_argument('-k',
                        '--kmer_len',
                        type=int,
                        action='store',
                        dest='k',
                        help='Specify the kmer/palindrome length.')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    This script receive a path to the csv files, a directory name to save the
    merged csv files and a list of species names.
    The final result is a merged csv file withh the mean of all kmer counts in the
    path directory.
    """
    print('Starting to process the script merge_csvs.py')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.path
    sub_dir = opt.sub_dir
    dir_out = opt.dir_out
    k = opt.k
    sub_sub_dir = f'kmer{k}'
    filenames = opt.filenames
    species = get_species_name(filenames)
    extension = opt.extension
    spc_data = defaultdict(list)
    for name in species:
        full_name_dir = os.path.join(dir_name, name, sub_dir, sub_sub_dir)
        spc_data[name] = spc_data.get(name,
                                      []) + [i for i in glob.glob(full_name_dir + f'/*.{extension}')]
    names_for_df = spc_data.keys()
    dfs = []
    cnt = 0
    for name in names_for_df:
        for filename in spc_data[name]:
            print(f'Working with the file: {filename}')
            dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
            df_final = reduce(lambda left, right: pd.merge(left, right, on='kmer'),
                              dfs).set_index('kmer').sort_index()
            name_csv = os.path.join(dir_out, name, sub_dir, name, sub_sub_dir)
            if os.path.exists(name_csv):
                pass
            else:
                make_me_a_folder(name_csv)
            df_final.to_csv(f'{name_csv}_all_concatenate.csv',
                            index=True)
            df_sum_rows = df_final.mean(axis=1)
            fullname = os.path.join(dir_out, name, sub_dir, name, sub_sub_dir)
            if os.path.exists(fullname):
                pass
            else:
                make_me_a_folder(fullname)
            df_sum_rows.to_csv(f'{fullname}_concatenate.csv',
                               index=True)
        cnt += 1
    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total files: {cnt}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!')


# df_merged = defaultdict(list)
# for name, df_list in dfs.items():
#     df_merged[name] = df_merged.get(name, [])
#     if len(df_list) > 1:
#         df_merged[name].append(reduce(lambda left, right: pd.merge(left,
#                                                                    right,
#                                                                    on='kmer'),
#                                       df_list).set_index('kmer').sort_index().mean(axis=1))
#     else:
#         df_merged[name].append(df_list)

# for name, data in df_merged.items():
#     name = name
#     df = pd.DataFrame(data).T
#     fullname = os.path.join(dir_name, name, subdir)
#     if not os.path.exists(fullname):
#         os.makedirs(fullname)
#     df.to_csv(f'{fullname}/{name}_kmer{k}_concatenate.csv', index=True)


if __name__ == "__main__":
    sys.exit(main())
