import sys
import os
import itertools
import multiprocessing
import time
from eggnogmapper import annota
import cPickle

def iter_tsv_lines(fname):
    for line in open(fname):
        if line.startswith('#') or not line:
            continue
        fields = map(str.strip, line.split('\t'))
        yield fields

def iter_hits(fname):
    last_query = None
    query_hits = []
    for fields in iter_tsv_lines(fname):
        if last_query is None:
            last_query = fields[0]
        elif fields[0] != last_query:
            yield last_query, query_hits
            last_query = fields[0]
            query_hits = []
        query_hits.append(fields)
    yield last_query, query_hits

def get_gos(members):
    all_gos = set()
    for m in members:
        all_gos.update(MEMBERGOS[m])
    all_gos.discard('')
    all_gos.discard('None')
    return all_gos

def annotate_blast_hits(arguments):
    t1 = time.time()
    annota.connect()
    query, query_hits = arguments

    all_hits = [(h[1], float(h[2]), float(h[10]), float(h[11])) for h in query_hits]
    blast_hits = set([str(h[0]) for h in all_hits if h[2] <= TARGET_CUTOFFS[0]])
    blast_hits.add(query)

    orthologs = QUERY_ORTHOLOGS[query]

    cutoff2result = {}
    for cutoff in TARGET_CUTOFFS:
        # get all hits that are not from self species and pass current cutoff
        filtered_blast_hits = set([str(h[0]) for h in all_hits
                                   if not h[0].startswith('%s.'%TARGET_TAXA)
                                   and h[2] <= cutoff])
        all_gos = get_gos(filtered_blast_hits)
        emapper_gos = get_gos(filtered_blast_hits & orthologs)
        cutoff2result[cutoff] = (query, all_gos, emapper_gos)

    return cutoff2result

def main(blast_hits_file):
    # Initialize output files
    OUT = {}
    for cutoff in TARGET_CUTOFFS:
        outfile = '%s.%s.%s.blast_annotations' %(blast_hits_file, GO_MODE, cutoff)
        emapper_outfile = 'emapper/%s/%s.%s.emapper.blast_filtered_annotations' %(EMAPPER_TAG, GO_MODE, cutoff)
        OUT[str(cutoff)] = [open(outfile, "w"), open(emapper_outfile, "w")]

    counter = 0
    for arguments in iter_hits(blast_hits_file):
        result = annotate_blast_hits(arguments)
        for cutoff, (query, all_gos, emapper_gos) in result.items():
            print >>OUT[str(cutoff)][0], '\t'.join([query, ','.join(all_gos)])
            print >>OUT[str(cutoff)][1], '\t'.join([query, ','.join(emapper_gos)])

        counter += 1
        if counter % 10 == 0:
            print "\r", counter,
            sys.stdout.flush()

    for f1, f2 in OUT.values():
        f1.close()
        f2.close()

if __name__ == "__main__":
    hits_file = sys.argv[1]
    HIT_IDS = set([n.strip() for n in open(sys.argv[2])])
    TARGET_TAXA = str(int(sys.argv[3]))
    GO_MODE = sys.argv[4]
    EMAPPER_TAG = sys.argv[5]
    CPU = int(sys.argv[6])

    sp, EMAPPER_METHOD, ORTHO_TYPE, TAX_SCOPE, GO_MODE, SELF = EMAPPER_TAG.split('.')
    EMAPPER_ORTHOLOGS_FILE = 'emapper/%s/%s.emapper.annotations.orthologs' %(TARGET_TAXA, EMAPPER_TAG)

    TARGET_CUTOFFS = [1E-03, 1E-10, 1E-40]
    GO_EXCLUDED =  set(["ND", "IEA"])
    if GO_MODE == 'experimental':
        GO_EVIDENCE = set(["EXP","IDA","IPI","IMP","IGI","IEP"])
    elif GO_MODE == "non-electronic":
        GO_EVIDENCE = None
    else:
        ValueError("Invalid go_mode")

    cache_file = "%s.%s.cache.pkl" %(hits_file, GO_MODE)
    try:
        print 'loading cache...'
        MEMBERGOS = cPickle.load(open(cache_file))
    except:
        print "Failed. Generating cache file"
        annota.connect()
        MEMBERGOS = annota.get_by_member_gos(HIT_IDS,
                                             target_go_ev=GO_EVIDENCE,
                                             excluded_go_ev=GO_EXCLUDED)
        cPickle.dump(MEMBERGOS, open(cache_file, "wb"), protocol=2)

    print len(MEMBERGOS), "hit ids in cache"

    main(hits_file)
