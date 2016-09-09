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
    false_positives = {}
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
        gos = set(map(str.strip, fields[6].split(',')))
        gos.discard('')
        emapper[query] = augment_gos(gos)
    return emapper

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

def benchmark(target_taxa, target_cutoffs, emapper_mode, emapper_tag):
    # Discard root GO ids (not really an annotation)
    root_levels = set(['GO:0008150', 'GO:0003674', 'GO:0005575'])

    false_positives = get_false_positives()
    # Check false positeves are loaded for all target taxa
    for i in TARGET_TAXA:
        false_positives[i]

    benchmark = {}
    for taxa in target_taxa:
        emapper_annot_file = "emapper_%s/%s/%s.%s.annot" %(emapper_mode, taxa, taxa, emapper_tag)
        interpro_annot_file = "interpro/%s.interpro.gz" %(taxa)
        print ' ', emapper_tag, emapper_annot_file
        print ' ', emapper_tag, interpro_annot_file
        emapper = get_augmented_emapper_gos(emapper_annot_file)
        interpro = get_augmented_interpro_gos(interpro_annot_file)

        print taxa, len(emapper), len(interpro)

        labels = []
        rows = []

        for cutoff in target_cutoffs:
            blast_annot_file = "emapper_%s/%s/blast_benchmark.%s.%s.annot.%s" %(emapper_mode, taxa, taxa, emapper_tag, cutoff)
            print ' ', emapper_tag, blast_annot_file
            matrix = []
            for line, fields in read_tsv(blast_annot_file):
                try:
                    (query, _, _, _, __htp, __otp,
                     tp_txt, hgos_txt, ogos_txt) = fields
                except:
                    continue

                tp = augment_gos(set(tp_txt.split(','))) - root_levels     # expected GOs (curated terms in query)
                hgos = augment_gos(set(hgos_txt.split(','))) - root_levels # blast GOs no selfhits (homologs)
                ogos = augment_gos(set(ogos_txt.split(','))) - root_levels # emapper GOS no selfhits (one2one orthologs)
                igos = interpro.get(query, set()) - root_levels # interpro results
                mgos = emapper.get(query, set()) - root_levels  # emapper results

                # sanitize sets with empty data
                for s in [tp, hgos, ogos, igos, mgos]:
                    s.discard('')

                htp = len(hgos & tp)
                otp = len(ogos & tp)
                itp = len(igos & tp)
                mtp = len(mgos & tp)

                # Sanity check: same tp numbers should have been inferred during processing.
                # You need to empty root_levels for this check to work!
                if not root_levels:
                    try:
                        assert htp == int(__htp)
                        assert otp == int(__otp)
                    except AssertionError:
                        print query, cutoff
                        print htp, __htp
                        print otp, __otp
                        raise

                # Solve conflicts between tp and fp definitions.
                fp = false_positives[taxa] - tp

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
    target_taxa, evalue_cutoffs, emapper_mode, emapper_tag = args
    #return emapper_tag, {emapper_tag:1}
    return emapper_tag, benchmark(target_taxa, evalue_cutoffs, emapper_mode, emapper_tag)

if __name__ == "__main__":
    TARGET_TAXA = ['9606', '7227', '3702', '4932', '511145']
    IDENT_CUTOFFS = ["20", "30", "40", "50", "60", "70"]
    EVALUE_CUTOFFS = ["0.001", "1e-10", "1e-40"]
    TARGET_TAXA_SCOPE = ["auto", "NOG"]
    TARGET_SELF = ['', 'excluded_self.']
    TARGET_ORTHO_TYPE = ['all', 'one2one']
    EMAPPER_MODES = ['hmm', 'diamond']
    BASEDIR = 'emapper'
    bench = {}
    cmds = []
    for emapper_mode in EMAPPER_MODES:
        for tself in TARGET_SELF:
            for tscope in TARGET_TAXA_SCOPE:
                for otype in TARGET_ORTHO_TYPE:
                    if emapper_mode == 'diamond':
                        emapper_tag = "dmnd.%s%s.%s" %(tself, otype, tscope)
                    else:
                        emapper_tag = "%s%s.%s" %(tself, otype, tscope)
                    cmds.append([TARGET_TAXA, EVALUE_CUTOFFS, emapper_mode, emapper_tag])

    pool = Pool(20)
    for tag, b in pool.imap(compute_benchmark, cmds):
        bench[tag] = b
        print tag, "Done"
        print b

    with open('all_benchmark_tables.pkl', 'w') as BENCH:
        cPickle.dump(bench, BENCH)
