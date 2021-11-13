#!/usr/bin/env python
# python csv_cleaning.py -di Results -sd kmer_count -do Results -k 12
import os
import sys
import time
import argparse
import csv
from collections import defaultdict


def get_csv_clean(csv_file):
    alpha = set('ACGT')
    km_cnt = defaultdict(float)
    with open(csv_file, 'r') as fh:
        data = csv.reader(fh)
        for row in data:
            kmer, cnt = row[0], float(row[1])
            if set(kmer).issubset(alpha):
                km_cnt[kmer] = km_cnt.get(kmer, 0.0) + cnt
    return km_cnt


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description="""A script to get the kmer frequency
     from csv files with kmer counts from genomes.""")

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Subdirectory name for output files.')  # kmer_count

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')  # Results/kmer_freq

    return parser.parse_args()


def main():
    """Parses options from the command line.
    Computes the k-mers frequencies from csv files with the kmer counts.
    """
    # checking the directory
    cwd = os.getcwd()
    print(f'The working directory: {cwd}')
    # counting time 
    start_time = time.process_time()
    # passing args
    arg = parse_arguments()
    sub_dir = arg.sub_dir
    dir_out = arg.dir_out
    file_amb = 'csv_to_clean'
    names_ambigous = defaultdict(str)
    with open(file_amb, 'r') as fh:
        for line in fh:
            name = line.strip().split('/')[2]
            names_ambigous[name] = names_ambigous.get(name, '')
            names_ambigous[name] += line.strip()
    print(f'number files: {len(names_ambigous)}')
    # checking if the output directory exist
    # if not make it
    f_pwd = os.path.join('Results', 'kmer_counts')
    # get the genus names
    cnt = 0
    for name, filename in names_ambigous.items():
        cleaned = get_csv_clean(filename)
        full_path = os.path.join(f_pwd, name)
        if os.path.exists(full_path):
            print(f'The path {full_path} exist')
            pass
        else:
            os.makedirs(full_path)
        csv_name = f'{full_path}/{name}_k2_8_chr.csv'
        print(f'Checking the full path {csv_name}')
        with open(csv_name, 'w') as fout:
            for km, cn in cleaned.items():
                fout.write(f'{km},{cn}\n')
        cnt += 1
    # get final time of the script
    end = time.process_time()
    total_time = end - start_time
    print(f'The script takes {total_time} to finish!')
    print(f'Where read and manipulated {cnt} files')
    print('Done!')


if __name__ == "__main__":
    sys.exit(main())
