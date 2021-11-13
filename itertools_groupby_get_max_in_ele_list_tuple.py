from itertools import groupby

groups = []
for k, g in groupby(sorted(data), key=lambda x: x[0]):
    groups.append(list(g))

for g in groups:
    print(g[0][0], 'min:', min(int(i[1]) for i in g), 'max:', max(int(i[2]) for i in g))
