from sys import argv

dur = None
try:
    dur = int(argv[1])
except:
    dur = 3 # default

from collections import defaultdict
graphs = defaultdict(dict)
measures = set()
with open(f'single_scalar_{dur}sec.dat') as data:
    for line in data:
        f = line.split()
        graph = f[-1] + f[0] + f[1] # replica generator order
        measure = f[2]
        value = None
        c = f[3]
        if c in ['True', 'False']:
            value = 1 * bool(c)
        elif 'j' in c:
            value = complex(c).real
        else:
            value = float(c)
        graphs[graph][measure] = value
        measures.add(measure)
M = list(measures)
print('graph', ' '.join(M))
for G in graphs:
    v = graphs[G]
    print(G, ' '.join([str(v.get(m, 'NA')) for m in M]))

        
    
