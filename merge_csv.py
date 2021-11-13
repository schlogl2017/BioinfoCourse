import os
import glob
import pandas as pd
import csv

os.chdir('\csv_files_direction')

extension = 'csv'
files = [i for i in glob.glob('*.{}'.format(extension))]
out_merg = ('\merged_csv_file_direction')
in_names = [pd.read_csv(f, delimiter=';', usecols = ['grid']) for f in files]
result = pd.concat(in_names)




