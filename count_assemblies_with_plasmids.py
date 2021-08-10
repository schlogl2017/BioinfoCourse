#!usr/bin/env python
import sys
from fasta_parser import fasta_item_counter, parse_fasta
from system_utils import get_fasta_files

if len(sys.argv) < 2:
    print('USAGE: < count_assemblies_with_plasmids.py > < directory name > ')
    sys.exit(1)


path = sys.argv[1]

assemblies_plasmids = []
cnt = 0
for filename in get_fasta_files(path):
    name = filename.split('/')[-1]
    headers = [header for header in parse_fasta(filename) if 'plasmid' in header]
    cnt += len(headers)
    assemblies_plasmids.append((set(headers), cnt))

print(f'The number of assemblies with plasmids are {cnt}')

with open('assemblies_with_plasmids.txt', 'w') as fo:
    for name in assemblies_plasmids:
        fo.write(f'{name}\n')
