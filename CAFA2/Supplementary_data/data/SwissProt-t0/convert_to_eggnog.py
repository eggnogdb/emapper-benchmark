import sys
from collections import defaultdict

all_gos = sys.argv[1]

name2gos = defaultdict(set)
for line in open(all_gos):
    name, go = line.strip().split()
    name2gos[name].add(go)

eggnog2gos = defaultdict(set)
inconsistent = set()
for line in open("uniprot2"):
    tax, seqname, sname, source = line.strip().split('\t')
    eggnog = "%s.%s" %(tax, seqname)
    if sname in name2gos:
        if eggnog in eggnog2gos and eggnog2gos[eggnog] != name2gos[sname]:
            #print eggnog2gos[eggnog] ^ name2gos[sname]
            #print eggnog, sname
            inconsistent.add(eggnog)
        else:
            eggnog2gos[eggnog] = name2gos[sname]

for k, v in eggnog2gos.iteritems():
    if k not in inconsistent:
        print '\t'.join([k, ','.join(sorted(v))])

    
