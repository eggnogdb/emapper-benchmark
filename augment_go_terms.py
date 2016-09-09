import os
import sys

BASE_DIR = os.path.split(os.path.abspath(__file__))[0]
GO2PARENTS = {}
for line in open('%s/go.tree' % (BASE_DIR)) :
    fields = map(str.strip, line.split('\t'))
    GO2PARENTS[fields[0]] = fields

def augment_gos(gos):
    augmented_gos = set()
    for g in gos:
        augmented_gos.update(GO2PARENTS.get(g, [g]))
    augmented_gos.discard('')
    return augmented_gos


if __name__ == "__main__":
    target_file = sys.argv[1]
    for line in open(target_file):
        fields = map(str.strip, line.split('\t'))
        augmented_gos = sorted(augment_gos(fields[1:]))
        print '\t'.join([fields[0], ','.join(augmented_gos)])
