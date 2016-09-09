import sys
import itertools
from eggnogmapper import annota

annota.connect()
CACHE = True

def iter_tsv_lines(fname):
    for line in open(fname):
        if line.startswith('#') or not line:
            continue
        fields = map(str.strip, line.split('\t'))
        yield fields

def main(blast_hits_file, emapper_hits_file, target_taxa, exec_mode, target_cutoffs):
    # pre-load emapper predictions
    seq2emapper = {}
    for fields in iter_tsv_lines(emapper_hits_file):
        #query_name || best_hit_eggNOG_ortholog || best_hit_evalue || best_hit_score
        #predicted_name (one-to-one) || orthologs (one-to-one) || GO (one-to-one)
        #KEGG_pathway (one-to-one) || predicted_name (all orthologs) || orthologs (not
        #one-to-one) || GO (not one-to-one) || KEGG_pathway (not one-to-one)

        qname = fields[0].strip()
        orthologs = set([o.strip() for o in fields[5].split(',')])
        orthologs.discard('')
        seq2emapper[qname] = orthologs

    # Initialize output files
    OUT = {}
    for cutoff in target_cutoffs:
        OUT[str(cutoff)] = open('blast_benchmark.%s.%s' %(mapper_file, cutoff), "w")

    query = None
    all_hits = []
    counter = 0
    missing_mapper_hits = 0

    # Process each block of blast hits. Assumes all hits for each query are
    # consecutive
    for query, query_hits in itertools.groupby(iter_tsv_lines(blast_hits_file), lambda x: x[0]):

        all_hits = [(h[1], float(h[2]), float(h[10]), float(h[11])) for h in query_hits]

        if query in seq2emapper:
            ortho_hits = seq2emapper[query]
        else:
            ortho_hits = set()
            missing_mapper_hits += 1

        blast_hits = set([str(h[0]) for h in all_hits])
        blast_hits.add(query)
        all_names = blast_hits | ortho_hits
        annotations = annota.get_by_member_annotations(all_names, set(["IEA", "ND"]))

        if CACHE:
            expected_gos = annotations[query][0]
        else:
            _, expected_gos, _ = annota.get_member_annotations([query], set(["IEA", "ND"]))

        for cutoff in target_cutoffs:
            # get all hits that are not from self species and pass current cutoff
            if exec_mode == "evalue":
                excluded_names = set([str(h[0]) for h in all_hits
                                        if h[0].startswith('%s.'%target_taxa)
                                        or h[2] > cutoff])

            elif exec_mode == 'identity':
                excluded_names = set([str(h[0]) for h in all_hits
                                    if h[0].startswith('%s.'%target_taxa)
                                    or h[1] < cutoff])

            filtered_blast_hits = blast_hits - excluded_names

            # sanity check: is there any predicted ortholog not part of the
            # blast result? Should not occur
            #if cutoff == 20:
            #    uncatched_orthologs = augmented_ortho_hits - blast_hits
            #    if uncatched_orthologs:
            #        pass
            #        #print len(uncatched_orthologs), "not seen", query, uncatched_orthologs

            filtered_ortho_hits = ortho_hits & filtered_blast_hits

            blast_gos = set()
            if filtered_blast_hits:
                if not CACHE:
                    _, blast_gos, _ = annota.get_member_annotations(filtered_blast_hits, set(["IEA", "ND"]))
                else:
                    [blast_gos.update(annotations[name][0]) for name in filtered_blast_hits]

            mapper_gos = set()
            if filtered_ortho_hits:
                if not CACHE:
                    _, mapper_gos, _ = annota.get_member_annotations(filtered_ortho_hits, set(["IEA", "ND"]))
                else:
                    [mapper_gos.update(annotations[name][0]) for name in filtered_ortho_hits]


            print >>OUT[str(cutoff)], '\t'.join(map(str, (query,
                                                          len(expected_gos),
                                                          len(blast_gos),
                                                          len(mapper_gos),
                                                          len(blast_gos & expected_gos),
                                                          len(mapper_gos & expected_gos),
                                                          ','.join(sorted(expected_gos)),
                                                          ','.join(blast_gos),
                                                          ','.join(mapper_gos),
                                                          )))

        counter += 1
        print counter, missing_mapper_hits


    for f in OUT.values():
        f.close()

if __name__ == "__main__":
    hits_file = sys.argv[1]
    mapper_file = sys.argv[2]
    target_taxa = sys.argv[3]
    exec_mode = sys.argv[4]

    if exec_mode == "identity":
        CUTOFFS = [20, 30, 40, 50, 60, 70]
    elif exec_mode == 'evalue':
        CUTOFFS = [1E-03, 1E-10, 1E-40]
    else:
        raise ValueError("choose mode: identity|evalue")

    main(hits_file, mapper_file, target_taxa, exec_mode, CUTOFFS)
