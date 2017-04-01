from collections import defaultdict
import cPickle
import gzip
from multiprocessing import Pool
import numpy as np
import pandas as pd
from augment_go_terms import augment_gos

def read_tsv(fname):
    if fname.endswith('.gz'):
        F = gzip.open(fname, 'r:gz')
    else:
        F = open(fname)

    for line in F:
        if line.startswith('#') or not line:
            continue
        yield line, map(str.strip, line.split('\t'))

def get_false_positives():
    false_positives = defaultdict(set)
    for line in open('data/all_negative_terms_per_taxid.tsv'):
        fields = map(str.strip, line.split('\t'))
        taxid = fields[0]
        if taxid == "83333":
            taxid = "511145"
        if taxid in false_positives:
            ValueError('Duplicate false positive entry %s' %fields[0])
        false_positives[taxid] = set(fields[1:])
    return false_positives

def get_augmented_emapper_gos(target_file):
    emapper = defaultdict(set)
    for line, fields in read_tsv(target_file):
        query = fields[0]
        gos = set(map(str.strip, fields[5].split(',')))
        gos.discard('')
        emapper[query] = augment_gos(gos)
    return emapper

def get_augmented_blast_gos(target_file):
    blast = defaultdict(set)
    for line, fields in read_tsv(target_file):
        query = fields[0]
        gos = set(map(str.strip, fields[1].split(',')))
        gos.discard('')
        blast[query] = augment_gos(gos)
    return blast

get_augmented_ref_gos = get_augmented_blast_gos

def get_augmented_interpro_gos(target_file):
    interpro = defaultdict(set)
    for line, fields in read_tsv(target_file):
        if "GO:" in line:
            query = fields[0]
            source = fields[3]
            gos = set(map(str.strip, fields[13].split('|')))
            gos.discard('')
            if not gos:
                raise ValueError('GOs in a different column?')
            interpro[query] = augment_gos(gos)
    return interpro

def benchmark(target_taxa, target_cutoffs, emapper_tag, blast_tag):
    # Discard root GO ids (not really an annotation)
    root_levels = set(['GO:0008150', 'GO:0003674', 'GO:0005575'])

    false_positives = get_false_positives()
    # Check false positeves are loaded for all target taxa
    for i in TARGET_TAXA:
        false_positives[i]

    benchmark = {}
    for taxa in target_taxa:
        interpro_annot_file = "interpro/%s.interpro.gz" %(taxa)
        if "non-electronic" in emapper_tag:
            groundtruth_file = "groundtruth/%s.groundtruth.non-electronic" %(taxa)
        elif "experimental" in emapper_tag:
            groundtruth_file = "groundtruth/%s.groundtruth.experimental" %(taxa)

        emapper_annot_file_self = "emapper/%s/%s.%s.self.emapper.annotations" %(taxa, taxa, emapper_tag)
        emapper_annot_file_noself = "emapper/%s/%s.%s.noself.emapper.annotations" %(taxa, taxa, emapper_tag)
        print ' ', emapper_tag, emapper_annot_file_noself
        print ' ', emapper_tag, emapper_annot_file_self
        print ' ', emapper_tag, interpro_annot_file

        emapper_self = get_augmented_emapper_gos(emapper_annot_file_self) # for interpro comparison
        interpro = get_augmented_interpro_gos(interpro_annot_file)
        groundtruth = get_augmented_ref_gos(groundtruth_file)

        print "Tax %s, eM(self) %s, iPro %s, truth %s" %(taxa, len(emapper_self), len(interpro), len(groundtruth))
        labels = []
        rows = []

        for cutoff in target_cutoffs:
            blast_annot_file = "blast/%s/%s.fa.hits.%s.%s.blast_annotations" %(taxa, taxa, blast_tag, cutoff)
            emapper_blast_file = "blast/%s/%s.%s.noself.%s.emapper.blast_filtered_annotations" %(taxa, taxa, emapper_tag, cutoff)

            blast = get_augmented_blast_gos(blast_annot_file)
            emapper_blast = get_augmented_blast_gos(emapper_blast_file)

            print ' ', emapper_tag, blast_annot_file
            matrix = []

            for query in blast:
                #try:
                #    (query, _, _, _, __htp, __otp,
                #     tp_txt, hgos_txt, ogos_txt) = fields
                #except:
                #    continue

                tp = groundtruth[query]
                hgos = blast[query]
                ogos = emapper_blast[query]
                igos = interpro[query]
                mgos = emapper_self[query]

                # Solve conflicts between tp and fp definitions.
                fp = false_positives[taxa] - tp

                #tp = augment_gos(set(tp_txt.split(','))) - root_levels     # expected GOs (curated terms in query)
                #hgos = augment_gos(set(hgos_txt.split(','))) - root_levels # blast GOs no selfhits (homologs)
                #ogos = augment_gos(set(ogos_txt.split(','))) - root_levels # emapper GOS no selfhits (one2one orthologs)
                #igos = interpro.get(query, set()) - root_levels # interpro results
                #mgos = emapper.get(query, set()) - root_levels  # emapper results

                # sanitize sets with empty data
                for s in [tp, fp, hgos, ogos, igos, mgos]:
                    s.discard('')

                htp = len(hgos & tp)
                otp = len(ogos & tp)
                itp = len(igos & tp)
                mtp = len(mgos & tp)

                hfp = len(hgos & fp)
                ofp = len(ogos & fp)
                ifp = len(igos & fp)
                mfp = len(mgos & fp)

                hunk = len(hgos - (fp | tp))
                ounk = len(ogos - (fp | tp))
                iunk = len(igos - (fp | tp))
                munk = len(mgos - (fp | tp))

                # TP/FP ratio over TP+FP (ignores unvalidated terms)
                htp_ratio = htp/float(htp+hfp) if htp+hfp else np.nan
                otp_ratio = otp/float(otp+ofp) if otp+ofp else np.nan
                itp_ratio = itp/float(itp+ifp) if itp+ifp else np.nan
                mtp_ratio = mtp/float(mtp+mfp) if mtp+mfp else np.nan

                matrix.append([htp_ratio, otp_ratio, itp_ratio, mtp_ratio,
                               htp,       otp,       itp,       mtp,
                               hfp,       ofp,       ifp,       mfp,
                               hunk,      ounk,      iunk,      munk,
                              ])

            # Converts matrix into a dataframe that we can use to operate
            header = map(str.strip, """
               htp_ratio    otp_ratio    itp_ratio   mtp_ratio
               htp          otp          itp         mtp
               hfp          ofp          ifp         mfp
               hunk         ounk         iunk        munk
               """.split())
            bench = pd.DataFrame(matrix, columns=header)

            # Summarize data and save stats line
            labels.append(cutoff)
            summary = [
                    bench.htp_ratio.mean(),
                    bench.otp_ratio.mean(),
                    bench.itp_ratio.mean(),
                    bench.mtp_ratio.mean(),

                    bench.htp.mean(),
                    bench.otp.mean(),
                    bench.itp.mean(),
                    bench.mtp.mean(),

                    bench.hfp.mean(),
                    bench.ofp.mean(),
                    bench.ifp.mean(),
                    bench.mfp.mean(),

                    bench.hunk.mean(),
                    bench.ounk.mean(),
                    bench.iunk.mean(),
                    bench.munk.mean(),

                    bench[(bench.htp > 0) & (bench.hfp + bench.hunk == 0)].htp.count(),
                    bench[(bench.otp > 0) & (bench.ofp + bench.ounk == 0)].otp.count(),
                    bench[(bench.itp > 0) & (bench.ifp + bench.iunk == 0)].itp.count(),
                    bench[(bench.mtp > 0) & (bench.mfp + bench.munk == 0)].mtp.count(),

                    bench[(bench.htp > 0) & (bench.hfp + bench.hunk > 0)].htp.count(),
                    bench[(bench.otp > 0) & (bench.ofp + bench.ounk > 0)].otp.count(),
                    bench[(bench.itp > 0) & (bench.ifp + bench.iunk > 0)].itp.count(),
                    bench[(bench.mtp > 0) & (bench.mfp + bench.munk > 0)].mtp.count(),

                    bench[(bench.htp == 0) & (bench.hfp + bench.hunk > 0)].htp.count(),
                    bench[(bench.otp == 0) & (bench.ofp + bench.ounk > 0)].otp.count(),
                    bench[(bench.itp == 0) & (bench.ifp + bench.iunk > 0)].itp.count(),
                    bench[(bench.mtp == 0) & (bench.mfp + bench.munk > 0)].mtp.count(),
                    ]

            rows.append(summary)


        header_main = map(str.strip, """
          htp_ratio   otp_ratio    itp_ratio   mtp_ratio
          htp         otp          itp         mtp
          hfp         ofp          ifp         mfp
          hunk        ounk         iunk        munk
          hp_tponly   op_tponly    ip_tponly   mp_tponly
          hp_tpplus   op_tpplus    ip_tpplus   mp_tpplus
          hp_notp     op_notp      ip_notp     mp_notp
        """.split())

        # Create summary benchmark data frame per taxa
        main_table = pd.DataFrame(rows, index=labels, columns=header_main)
        benchmark[taxa] = main_table
        print taxa, "finished"

    return benchmark

def compute_benchmark(args):
    target_taxa, evalue_cutoffs, emapper_tag, blast_tag = args
    return emapper_tag, benchmark(target_taxa, evalue_cutoffs, emapper_tag, blast_tag)

if __name__ == "__main__":
    TARGET_TAXA = ['9606', '7227', '3702', '4932', '511145', '5833', '759272']
    EVALUE_CUTOFFS = ["0.001", "1e-10", "1e-40"]
    TARGET_TAXA_SCOPE = ["auto", "NOG"]
    TARGET_ORTHO_TYPE = ['all', 'one2one']
    EMAPPER_MODES = ['hmmer', 'diamond']
    TARGET_GO = ['non-electronic', 'experimental']

    bench = {}
    cmds = []
    for emapper_mode in EMAPPER_MODES:
        for tscope in TARGET_TAXA_SCOPE:
            for otype in TARGET_ORTHO_TYPE:
                for gtype in TARGET_GO:
                    emapper_tag = "%s.%s.%s.%s" %(emapper_mode, otype, tscope, gtype)
                    blast_tag = "%s" %(gtype)
                    cmds.append([TARGET_TAXA, EVALUE_CUTOFFS, emapper_tag, blast_tag])

    pool = Pool(20)
    for tag, b in pool.imap(compute_benchmark, cmds):
        bench[tag] = b
        print tag, "Done"
        print b

    with open('all_benchmark_tables_new.pkl', 'w') as BENCH:
        cPickle.dump(bench, BENCH)
