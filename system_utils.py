#!usr/bin/env python
# -*- coding: utf-8 -*-


import os
import glob
from collections import defaultdict, Counter
import collections
import pandas as pd
import dask.dataframe as dd
import numpy as np


def get_files(dir_name):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_names(path):
    return [(file.split('/')[3][:-3]) for file in get_files(path)]


def get_dir_names(filenames):
    dirn = set()
    subd = set()
    subsub = set()
    for filename in filenames:
        dir_splt = filename.split('/')
        dirn.add(dir_splt[0])
        subd.add(dir_splt[1])
        subsub.add(dir_splt[2])
    return next(iter(dirn)), next(iter(subd)), next(iter(subsub))


def get_dir_name(filename):
    subd, subsub = set(), set()
    names = filename.split('/')
    for i in range(len(names)):
        subd.add(names[1])
        subsub.add(names[2])
    return next(iter(subd)), next(iter(subsub))


def make_me_a_folder(folder_name):
    os.getcwd()
    try:
        os.makedirs(folder_name)
    except OSError:
        print("Directory already exists")
    return os.path.exists(folder_name)


def get_full_name(dir_root, sud_dir, ssub):
    return os.path.join(dir_root, sud_dir, ssub)


def get_fasta_files(dir_name):
    # tuple with the file extention to check
    ext = tuple(['fa.gz', '.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    infiles = []
    # lokking for file in the path
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            # getting the filenames
            input_files = os.path.join(path, name)
            if input_files.endswith(ext):
                infiles.append(input_files)
    return infiles


def get_files2(dir_name, subdir):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry, subdir)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_list_paths(spc_names, path, ssub):
    pwd = defaultdict()
    for name in spc_names:
        pwd[name] = get_full_name(path, name, ssub)
    return pwd


def get_fasta_files_paths(dict_paths):
    dict_fasta_paths = defaultdict(list)
    for name, path in dict_paths.items():
        dict_fasta_paths[name] = dict_fasta_paths.get(name, []) + glob.glob(path)
    return dict_fasta_paths


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def get_fasta_files_dict(dir_name):
    ext = tuple(['fa.gz', '.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    infiles = defaultdict(list)
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            infiles[name] = infiles.get(name, [])
            input_files = os.path.join(path, name)
            if input_files.endswith(ext):
                infiles[name].append(input_files)
    return infiles


def fasta_file_paths(fasta_dirs):
    fasta_paths = defaultdict(list)
    exts = tuple(['fa.gz', '.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    for name, pwd in fasta_dirs.items():
        fasta_paths[name] = fasta_paths.get(name, [])
        for path, _, files in os.walk(pwd):
            for n in files:
                input_files = os.path.join(path, n)
                if input_files.endswith(exts):
                    fasta_paths[name].append(input_files)
    return fasta_paths


def get_names_list(text_file):
    """
    Function to get the names of species or genus in a text file.
    
    Inputs:
    
        text_file - a string representing a text file where the 
                    species or genus are listed one by line.
    
    Outputs:
    
        names_list - a list/array-like that countain the name of the
                     species or genus.
                
        > get_names_list('teste_genus.txt')
        
        ['Acidiphilium']
    """
    # initialize the list
    names_list = []
    # open the text file
    with open(text_file, 'r') as fh:
        # iterates through the file handle
        for line in fh:
            # clean the spaces
            name = line.strip()
            # add the name to the list
            names_list.append(name)
    return names_list


def find_csv_filenames(dir_name, spc_name, sub_dir, ssub, num_k):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, f'{ssub}{str(i)}') for i in range(1, num_k)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def get_csvs_to_df_concatenation(spc_name, filenames_dict):
    """Function tha receives a species name and a dictionary with
    all csvs files to be concatenated.
    """
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename, compression='gzip', dtype={'kmer': str, 'count': np.int32}))
    dfs = [df.set_index("kmer", drop=True) for df in dfs]
    concat = pd.concat(dfs, axis=0, copy=False).reset_index()
    return concat, len(dfs)


def csv_filenames(dir_name, spc_name, sub_dir, ssub):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, ssub)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def path_generator(dir_name, sub_dir, genus, sub_sub_dir, type_file):
    """
    Function to generate one by one the complete paths to the 
    fasta files.
    The directory need to be designed like this:
        base directory/sub_directory/genus/sub_sub_directory/files.type
        Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz
        
    Inputs:
        
        dir_name - string representing the base directory.
        sub_dir - string representing the sub directory.
        genus - string representing the bacterial genus.
        sub_sub_dir - string representing the sub sub directory, where the files
                      are located.
        type_file - represents the final file extention (csv, gz, etc)
    
    Outputs:
    
        Yields the complete path to the source files.
        
        >list(path_generator('Data', 'Genomes_splitted', 'Acidiphilium', 'Chromosomes', 'gz'))
        
        ['Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz',
        'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz']
        
    """
    # generates the paths to the directory
    pwd = os.path.join(dir_name, sub_dir, genus, sub_sub_dir)
    # iterates through the directory
    for filename in os.listdir(pwd):
        # checks if the files have the correct type
        # and if so yields the paths
        if filename.endswith(type_file):
            yield os.path.join(pwd, filename)


def data_generator(file_path_dict):
    """
    Yields a tuples as names/genus and the full path to the source
    files.
    
    Inputs:
    
        file_path_dict - a dictionary-like object mapping the names/genus to
                         the source of files to be analized.
    
    Outputs:
    
        Yields a tuple with first element as the name or genus and the second
        element the full path to the source file.
    
    for data in data_generator(file_dict):
        print(data)
    
    ('Acidiphilium', 'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz')
    ('Acidiphilium', 'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz')
    """
    # iterates through the dictionary items
    for name, filenames in file_path_dict.items():
        # gets the name
        spc = name
        # and then iterates through the filenames
        # that is a list with the full path
        # then yields the name and the full path
        for filename in filenames:
            yield spc, filename


def file_dict(text_file, dir_name, sub_dir, sub_sub_dir, type_file):
    """
    Function to generate a dictionary-like object with the name/genus to
    the complete path to the fasta files.
    The directory need to be designed like this:
        base directory/sub_directory/genus/sub_sub_directory/files.type
        Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz
        
    Inputs:
    
        text_file - a string representing a text file where the 
                    species or genus are listed one by line.        
        dir_name - string representing the base directory.
        sub_dir - string representing the sub directory.
        sub_sub_dir - string representing the sub sub directory, where the files
                      are located.
        type_file - represents the final file extention (csv, gz, etc)
    
    Outputs:
    
        paths_to_files - a dictionary-like object where keys are names and
                         values are the full path to the files

    >file_dict('teste_genus.txt', 'Data', 'Genomes_splitted', 'Chromosomes', 'gz')
    defaultdict(list,
            {'Acidiphilium': 
            ['Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000016725.1_chr.fna.gz',
              'Data/Genomes_splitted/Acidiphilium/Chromosomes/GCF_000202835.1_chr.fna.gz']})
    """
    # initialze the dictionary with default values as lists
    paths_to_files = defaultdict(list)
    # get the names to the keys
    names = get_names_list(text_file)
    # iterates through the names list
    for name in names:
        # add the keys to the dictionary
        # add the full paths to the values items
        # as a list and return the dictionary
        paths_to_files[name] = paths_to_files.get(name, []) + list(path_generator(dir_name,
                                                                                  sub_dir,
                                                                                  name,
                                                                                  sub_sub_dir,
                                                                                  type_file))
    return paths_to_files


def get_species_name(filenames):
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def get_csv_files(dir_name):
    # tuple with the file extention to check
    ext = tuple(['csv.gz', '.csv'])
    infiles = []
    # lokking for file in the path
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            # getting the filenames
            input_files = os.path.join(path, name)
            if input_files.endswith(ext):
                infiles.append(input_files)
    return infiles


def get_csvs_to_df_merge(spc_name, filenames_dict):
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
    dfs = [df.set_index("kmer", drop=True) for df in dfs]
    merged = pd.concat(dfs,
                       axis=1,
                       keys=range(len(dfs)),
                       copy=False).sort_index().mean(axis=1).astype(int).reset_index()
    return merged.rename(columns={'index': 'kmer', 0: 'count'}), len(dfs)


def dask_csv_merge(name, filenames_dict):
    """
    Function to read and merge a bunch of large csv files.
    Each csvs has two coluns 'kmer, count' and the final
    result is a pandas.core.series.Series. This object can be
    transformed in a pandas daraframe or a csv file.

    Inputs:
        name - a string representing a dictionary key.
        filename_dict - a dicitionary-like object mapping
                        the genus name to the csv file locations(as a list).

    Outputs:
        dfm - a pandas.core.series.Series where kmers (strings length k) are
              the indexes and the mean of all csv counts are the values.
    """
    dfs = [dd.read_csv(filename,
                       dtype={'kmer': str, 'count': np.int32},
                       blocksize=None).set_index('kmer') for filename in filenames_dict[name]]
    dfm = dd.concat(dfs, axis=1).mean(axis=1).astype(np.int32).compute()
    return dfm, len(dfs)


def concatenation_dfs(spc_name, dict_csvs_paths):
    csvs = dict_csvs_paths[spc_name]
    concat_df = pd.DataFrame()
    for csv in csvs:
        temp_df = pd.read_csv(csv, dtype={'kmers': str,
                                          'count': np.int32}).set_index("kmers",
                                                                        drop=True)
        concat_df = pd.DataFrame(pd.concat([concat_df,
                                            temp_df],
                                           axis=1,
                                           join='outer',
                                           copy=False).sort_index().mean(axis=1))
        del temp_df
    concat_df[0] = concat_df[0].astype('float32')
    return concat_df.reset_index().rename(columns={'index': 'kmers',
                                                   0: 'counts'}), len(csvs)


def list_to_string(org_list, seperator=' '):
    """ Convert list to string, by joining all item in list with given separator.
        Returns the concatenated string """
    return seperator.join(org_list)


def seq_to_dict(dict_paths):
    """
    This function returns a dictionary-like object with sequence ids as keys and
    sequence (as strings) as values.
    
    Inputs:
        
        dict_paths - a dictionary-like object with bacterial genus as keys and 
                     a list of fasta file paths as values.
        
    Outputs:
    
        A tuple like object:
         
             seq_dict - a dictionary-like object where keys are sequence ids and
                        values are the sequence itself.
                        
             gen - bacterial genus name.
    
             total_seqs - number of the sequences that where found for that genus.
    
    """
    gen_name = ''
    total_seqs = 0
    dict_fasta = defaultdict(list)
    for gen, paths in dict_paths.items():
        gen_name += gen
        for path in paths:
            for seq_id, seq in parse_fasta(path):
                seq_id = '_'.join(name.split(' ')[:3])
                dict_fasta[seq_id] += [seq]
            total_seqs += 1
    return {k: list_to_string(v, '') for k, v in dict_fasta.items()}, gen_name, total_seqs


def get_sequences_from_dict(dict_paths):
    dict_fasta = defaultdict(list)
    for gen, paths in dict_paths.items():
        c = 0
        for path in paths:
            for rec in SeqIO.parse(gzip.open(path, 'rt'), 'fasta'):
                dict_fasta[rec.id] += [str(rec.seq)]
            c += 1
        print('Number seq:', c)
    return {k: list_to_string(v[0], '') for k, v in dict_fasta.items()}


# def get_sequence_lengths(fasta_filenames):
#     """Returns dictionary of sequence lengths, keyed by organism.
#     Biopython's SeqIO module is used to parse all sequences in the FASTA
#     file corresponding to each organism, and the total base count in each
#     is obtained.
#     NOTE: ambiguity symbols are not discounted.
#     """
#     tot_lengths = {}
#     for fn in fastafilenames:
#         tt_lengths[os.path.splitext(os.path.split(fn)[-1])[0]] = \
#             sum(o[len(s) for s in SeqIO.parse(fn, 'fasta')])
#     return tot_lengths


def save_csv_from_dict(kmer_counts, dir_out, spc, name, sub_sub_dir, sub_sub_sub_dir, max_k):
    df = pd.DataFrame(kmer_counts.items(), columns=['kmer', 'count'])
    full_name_dir = os.path.join(dir_out, spc, sub_sub_dir, sub_sub_sub_dir)
    csv_name = os.path.join(full_name_dir, name)
    if not os.path.exists(full_name_dir):
        os.makedirs(full_name_dir)
    print(f'The final csv file was saved in {full_name_dir}\n')
    df.to_csv(f'{csv_name}_{sub_sub_sub_dir}_{max_k}.csv', index=False)


def get_files_paths(dir_name, sub_dir_name):
    """Returns a dictionary like object using baterial genera
    as keys and a list of path to fasta fileas as values.
    Inputs:
        dir_name: directory name
        sub_dir_name : sub directory name
    Outputs:
        dictionary: with genera name and a list of 
                    files.
    
    Example:
    dir is like: Data/bacteria_splitted/Mycolicibacillus/chromosomes
                 Data/bacteria_splitted/Mycolicibacillus/plasmids (if the case)
    dirname: 'Data/bacteria_splitted'
    sub_dir: 'chromosomes'
    fasta_dicts = get_files(dir_name, 'chromosomes')
    fasta_dicts['Mycolicibacillus']
    ['NZ_AP022594.1_Mycolicibacillus_koreensis_strain_JCM_19956.fna.gz']
    """
    # create a list of file and sub directories
    # names in the given directory
    spc_names = sorted(os.listdir(dir_name))
    all_files = defaultdict(list)
    for name in spc_names:
        all_files[name] = all_files.get(name, [])
        paths = os.path.join(dir_name, name, sub_dir_name)
        for file in os.listdir(paths):
            full_paths = os.path.join(paths, file)
            all_files[name].append(full_paths)
    return all_files


def get_name_path(path_dict):
    for name, pwd in path_dict.items():
        yield name, pwd


def get_all_fasta(dir_in, ssubdir_in, extension):
    """
    Function to fecth all fasta files from a directory.
    The directory has a tree like:
    root/subdir/genus_name/subsub_dir. Ex: 'Data/Genomes_splitted/Salmonella/Chromosomes'.

    Inputs:
        dir_in - a string representing the root directory.
        name - a string representing representing the name of a subdirectory,
               that represent a species or genus.
        ssubdir_in - a string representing a subdirectory inside the root directory.
        extension - a string representing the type of the file (compressed or not). Ex., 'gz'.

    Outputs:
        all_fasta -  a dictionary-like object mapping the genus name (as key) to their
                    list of full path to the fasta files.
    """
    names = sorted(os.listdir(dir_in))
    # initialize the container
    all_fasta = defaultdict(list)
    # get the full paths to the fasta files
    for name in names:
        full_path = os.path.join(dir_in, name, ssubdir_in)
        # add the names]keys and the fullpaths to the container
        all_fasta[name] = all_fasta.get(name, []) + glob.glob(full_path + f'/*.{extension}')
    return all_fasta
