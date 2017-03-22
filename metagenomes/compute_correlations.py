from collections import Counter
from eggnogmapper import annota
import re
import numpy as np
from scipy import stats

def get_gos(gene):
    go_excluded = set(["ND", "IEA"])
    pnames, gos, kegg = annota.get_member_annotations([gene], excluded_go_ev=None, target_go_ev=None)
    return gos

    
def get_stats(read_gos, emapper_gos, ipro_gos):
    govectors = {}
    govectors['GO terms in either prediction or truth'] = [set(read_gos.keys()) | set(emapper_gos.keys()),
                             set(read_gos.keys()) | set(ipro_gos.keys())]
    govectors['GO terms in both prediction or truth'] = [set(read_gos.keys()) & set(emapper_gos.keys()),
                              set(read_gos.keys()) & set(ipro_gos.keys())]
    truth = set(read_gos.keys())
    
    for name, vector in govectors.items():
        print name
        refe = []
        emapper = []
        for g in vector[0]:
            emapper.append(emapper_gos[g])
            refe.append(read_gos[g])
        print " Spearman: groundtruth vs emapper  R=%g e-value=%g" %(stats.spearmanr(refe, emapper))

        refi = []
        ipro = []
        for g in vector[1]:
            ipro.append(ipro_gos[g])
            refi.append(read_gos[g])
        print " Spearman: groundtruth vs interpro R=%g e-value=%g" %(stats.spearmanr(refi, ipro))

    ref = []
    emapper = []
    ipro = []
    for g in truth:
        emapper.append(emapper_gos[g])
        ipro.append(ipro_gos[g])
        ref.append(read_gos[g])

    print "GO terms only in groundtruth (%d):" %(len(truth))
    print " Spearman: groundtruth vs emapper  R=%g e-value=%g" %(stats.spearmanr(ref, emapper))
    print " Spearman: groundtruth vs interpro R=%g e-value=%g" %(stats.spearmanr(ref, ipro))


annota.connect()
for sample in range(4):
    read_gos = Counter()
    ipro_gos = Counter()
    emapper_gos = Counter()

    # Load abundances of predicted genes in the catalog
    gene2count = {}
    for line in open("simulated_samples/sample-%d.predicted_gene_abundances" %sample):
        if not line.strip() or line.startswith('#') or line.startswith("-1") or line.startswith("gene"):
            continue
        gene, count = map(str.strip, line.split('\t'))
        gene2count[gene] = float(count)

    # Load emapper annotations for predicted genes in the catalog
    for line in open("emapper_annotations/sample-%d.catalog.faa.emapper.annotations" %sample):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split('\t')
        qname = fields[0]
        gos = set([g.strip() for g in fields[5].split(',') if g.startswith('GO:')])
        for g in gos:
            emapper_gos[g] += gene2count[qname]

    # Load ipro annotations for predicted genes in the catalog
    for line in open("interpro_annotations/sample-%d.catalog.ipro" %sample):
        if not line.strip() or line.startswith('#'):
            continue
        if 'GO:' not in line:
            continue
        fields = line.split('\t')
        qname = fields[0]
        gos = [g.strip() for g in fields[13].split('|') if g.startswith('GO:')]
        for g in gos:
            ipro_gos[g] += gene2count[qname]

    # Load groundtruth annotations based on real gene abundances in the simulated set
    for line in open("simulated_samples/sample-%d.original_gene_abundances" %sample):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split('\t')
        count = float(fields[1])
        gos = get_gos(fields[0])
        for g in gos:
            read_gos[g] += count

    print "### Sample-%s GOs from reads:% 5d   GOs from emapper:% 5d   GOs from ipro:% 5d" %(sample, len(read_gos), len(emapper_gos), len(ipro_gos))
    print '```'
    get_stats(read_gos, emapper_gos, ipro_gos)
    print '-'*100
    print '```'
