#!usr/bin/env python
# Finding and counting the frequency of known pairs of words in multiple files - Stack Overflow


import os
import re
from collections import Counter


def pairs(text):
    ans = re.findall(r'[A-Za-z]+', text)
    return (tuple(ans[i:i+2]) for i in range(len(ans)-1))


mypairs = tuple([ tuple(line.split()[-2:]) for line in open('results.txt')])

c = Counter()
folderpath = 'path/to/directory'
for dirpath, dnames, fnames in os.walk(folderpath):
    for f in fnames:
        if not '.txt' in f:
            continue
        for line in open(os.path.join(dirpath, f)):
            c += Counter(p for p in pairs(line) if p in mypairs)
for item in c.most_common():
    print(item)
