import os
from matplotlib import pyplot as plt
from StringIO import StringIO
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
import ipywidgets as widgets
from IPython.display import display, Markdown, Latex

# CONFIG
PROTEOME_SIZE = {'9606': 22834, '511145': 4146, '7227': 13937, '3702': 28128, '4932':6692, '5833':5429, '759272':7402}
TARGET_TAXA = ['9606', '7227', '3702', '4932', '511145']
RESULTS_DIR = './'
EVALUE_CUTOFFS = ["0.001", "1e-10", "1e-40"]
ncbi = NCBITaxa()
TAXA_NAMES = {k:' '.join(v.split()[:2]) for k,v in ncbi.get_taxid_translator(PROTEOME_SIZE.keys()).items()}

pagebreak = "\\clearpage"

def plot_blast_benchmark(benchmark, target_taxa, target_cutoffs, ylabel, emapper_tag, nodisplay=False):
    f, axes = plt.subplots(5, 3, figsize=(16.9, 15))
    for i, taxa in enumerate(target_taxa):
        if taxa not in benchmark:
            continue
        data = benchmark[taxa]
        first_row = data.loc[target_cutoffs[0]]

        ylim = (-2, (len(target_cutoffs))*len(target_taxa) )
        ytick_labels = reversed(target_cutoffs)
        yticks = range(0, len(target_cutoffs)*len(target_taxa), 5)

        # PLOT FIRST COLUMN
        # =========================
        drate = axes[i][0]
        drate.set_axis_bgcolor('white')
        #drate.get_yaxis().set_visible(False)
        drate.spines['right'].set_visible(False)
        drate.spines['top'].set_visible(False)
        drate.spines['left'].set_visible(False)
        drate.get_xaxis().tick_bottom()
        drate.get_yaxis().tick_left()
        drate.set_xlabel("% TP and FP GO terms per protein", {'fontsize': 10})
        drate.set_yticks(yticks)
        drate.set_yticklabels(ytick_labels)
        drate.set_xlim((0, 100))
        drate.set_ylim(ylim)
        #drate.yaxis.tick_right()

        rn = 0
        for name, row in reversed([r for r in data.iterrows()]):
            htp_prop = row.htp_ratio * 100
            hfp_prop = 100 - htp_prop

            otp_prop = row.otp_ratio * 100
            ofp_prop = 100 - otp_prop


            a = drate.barh((rn-2), htp_prop, height=1.9, color='darkgreen',linewidth=0, alpha=0.4,
                           label="True Positives (Blast)")
            b = drate.barh((rn-2), hfp_prop, height=1.9, left=htp_prop,
                           color='indianred',linewidth=0, alpha=0.4,
                           label="False Positives (Blast)")

            c = drate.barh((rn+0.1), otp_prop, height=1.9, color='darkgreen', linewidth=0, alpha=0.9,
                           label="True Positives (eggNOG)")
            d = drate.barh((rn+0.1), ofp_prop, height=1.9, left=otp_prop,
                           color='indianred',linewidth=0, alpha=0.9,
                           label="False Positives (eggNOG)")

            rn += 5.0

        # PLOT SECOND COLUMN
        # =========================

        cov = axes[i][1]
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()

        cov.set_xlabel("Average GO terms per protein", {'fontsize': 10})
        cov.set_yticks([])
        cov.set_ylim(ylim)

        rn = 0
        max_gos = first_row.htp + first_row.hfp + first_row.hunk
        cov.set_xlim(0, max_gos)
        for name, row in reversed([r for r in data.iterrows()]):
            La = cov.barh((rn+0.1), row.otp,
                         height=1.9, color='darkgreen',linewidth=0, alpha=0.9,
                         label='True Positive terms (eggNOG)')
            Lb = cov.barh((rn+0.1), row.ofp, left=row.otp,
                         height=1.9, color='indianred',linewidth=0, alpha=0.9,
                         label='False Positive terms (eggNOG)')
            Lc = cov.barh((rn+0.1), row.ounk, left=(row.otp+row.ofp),
                         height=1.9, color='darkgrey',linewidth=0, alpha=0.9,
                         label='Not curated terms (eggNOG)')

            Ld = cov.barh((rn-2), row.htp,
                         height=1.9, color='darkgreen',linewidth=0, alpha=0.4, label='(BLAST)')
            Le = cov.barh((rn-2), row.hfp, left=row.htp,
                         height=1.9, color='indianred',linewidth=0, alpha=0.4, label='(BLAST)')
            Lf = cov.barh((rn-2), row.hunk, left=(row.htp+row.hfp),
                         height=1.9, color='darkgrey',linewidth=0, alpha=0.4, label='(BLAST)')


            hratio = row.htp / (row.htp + row.hfp + row.hunk)
            oratio = row.otp / (row.otp + row.ofp + row.ounk)
            cov.text(max_gos, rn+0.4, '%0.2f ' %oratio, horizontalalignment='right', color = "blue")
            cov.text(max_gos, rn-1.4, '%0.2f ' %hratio, horizontalalignment='right')
            rn += 5.0

            if taxa == '9606' and name == target_cutoffs[0]:
                cov.text(max_gos, rn-2.5, 'TP ratio\n(higher is better)',
                         horizontalalignment='right', fontweight="bold", fontsize=8)


        # PLOT THIRD COLUMN
        # =========================
        cov = axes[i][2]
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()
        cov.set_xlabel("% of proteome annotated", {'fontsize': 10})
        cov.set_yticks([])
        cov.set_ylim(ylim)
        #cov.get_xaxis().tick_right()

        rn = 0
        psize = PROTEOME_SIZE[taxa]/100
        cov.set_xlim(0, 100)
        for name, row in reversed([r for r in data.iterrows()]):
            Lg = cov.barh((rn+0.1), row.op_tponly/psize,
                         height=1.9, color='darkblue',linewidth=0, alpha=0.9,
                         label='Proteins with only TP terms (eggNOG)')

            Lh = cov.barh((rn+0.1), row.op_tpplus/psize, left=row.op_tponly/psize,
                         height=1.9, color='purple',linewidth=0, alpha=0.9,
                         label='Proteins with TP and other terms (eggNOG)')

            Li = cov.barh((rn+0.1), row.op_notp/psize, left=(row.op_tponly + row.op_tpplus)/psize,
                         height=1.9, color='orange',linewidth=0, alpha=0.9,
                         label='Proteins lacking TP terms (eggNOG)')

            Lj = cov.barh((rn-2), row.hp_tponly/psize,
                         height=1.9, color='darkblue',linewidth=0, alpha=0.4, label='(BLAST)')
            Lk = cov.barh((rn-2), row.hp_tpplus/psize, left=row.hp_tponly/psize,
                         height=1.9, color='purple',linewidth=0, alpha=0.4, label='(BLAST)')
            Lm = cov.barh((rn-2), row.hp_notp/psize, left=(row.hp_tponly + row.hp_tpplus)/psize,
                         height=1.9, color='orange',linewidth=0, alpha=0.4, label='(BLAST)')

            hratio = row.hp_notp / (row.hp_notp + row.hp_tponly + row.hp_tpplus)
            oratio = row.op_notp / (row.op_notp + row.op_tponly + row.op_tpplus)

            cov.text(100, rn+0.4, '%0.2f ' %oratio, horizontalalignment='right', color = "blue")
            cov.text(100, rn-1.4, '%0.2f ' %hratio, horizontalalignment='right')

            rn += 5.0
            if taxa == '9606' and name == target_cutoffs[0]:
                cov.text(100, rn-2.5, 'no-TP ratio\n(lower is better)',
                         horizontalalignment='right', fontweight="bold", fontsize=8)


        if i == 0:
            drate.legend(handles=[La, Lb, Lc, Ld, Le, Lf], loc=2, ncol=2, bbox_to_anchor=(0.1, 1.4), fontsize=10)
            cov.legend(handles=[Lg, Lh, Li, Lj, Lk, Lm], loc=2, ncol=2, bbox_to_anchor=(-0.9, 1.4), fontsize=10)
        #cov.yaxis.tick_right()
        cov.yaxis.set_label_position("right")
        drate.set_ylabel(ylabel, x=2, fontsize=10)
        cov.set_ylabel(TAXA_NAMES[int(taxa)], labelpad=30, fontsize=18, rotation=-90, style='italic')
        cov.get_yaxis().tick_right()
        # cov.margins(left=50)


    plt.subplots_adjust(left=None, bottom=0, right=None, top=1.02, wspace=None, hspace=None)

    plt.savefig("plots/emapper_vs_blast.%s.pdf" %(emapper_tag), facecolor='w', edgecolor='w',
            orientation='portrait', bbox_inches = 'tight')

    plt.savefig("plots/emapper_vs_blast.%s.png" %(emapper_tag), dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', bbox_inches = 'tight')
    if not nodisplay:
        plt.show()


def plot_interpro_benchmark(benchmark, target_taxa, emapper_tag, nodisplay=False):
    f, axes = plt.subplots(1, 3, figsize=(16.9, 5))
    max_gos = max([(benchmark[t].mtp + benchmark[t].mfp + benchmark[t].munk)[0] for t in target_taxa])

    rn = 0
    labels = []
    for taxa in reversed(target_taxa):
        labels.append(' '.join(TAXA_NAMES.get(int(taxa), taxa).split()[:2]))

        if taxa not in benchmark:
            rn += 5
            continue

        row = benchmark[taxa].loc["0.001"]
        # PLOT FIRST COLUMN
        # =========================
        drate = axes[0]
        yticks = range(2, len(target_taxa)*5, 5)
        ylim = (0, len(target_taxa)*5)
        if taxa == "9606":
            #drate.get_yaxis().set_visible(False)
            drate.spines['right'].set_visible(False)
            drate.spines['top'].set_visible(False)
            drate.spines['left'].set_visible(False)
            drate.get_xaxis().tick_bottom()
            drate.get_yaxis().tick_left()
            drate.set_xlabel("% TP and FP GO terms per protein", {'fontsize': 10})
            drate.set_yticks(yticks) #[2, 7, 12, 17, 22])
            drate.set_yticklabels(labels, fontsize=15, style='italic')
            drate.set_xlim((0, 100))
            drate.set_ylim(ylim)
            #drate.yaxis.tick_right()



        itp_prop = row.itp_ratio * 100
        ifp_prop = 100 - itp_prop

        otp_prop = row.mtp_ratio * 100
        ofp_prop = 100 - otp_prop


        c = drate.barh((rn+2), otp_prop, height=1.9, color='darkgreen', linewidth=0, alpha=0.9,
                       label="True Positives (eggNOG)")
        d = drate.barh((rn+2), ofp_prop, height=1.9, left=otp_prop,
                       color='indianred',linewidth=0, alpha=0.9,
                       label="False Positives (eggNOG)")

        a = drate.barh((rn-0.1), itp_prop, height=1.9, color='darkgreen',linewidth=0, alpha=0.4,
                       label="True Positives (Blast)")
        b = drate.barh((rn-0.1), ifp_prop, height=1.9, left=itp_prop,
                       color='indianred',linewidth=0, alpha=0.4,
                       label="False Positives (Blast)")


        # ROW IN SECOND COLUMN
        # =========================

        cov = axes[1]
        cov.set_xlabel("Average GO terms per protein", {'fontsize': 10})
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()
        cov.set_yticks([])
        cov.set_ylim(ylim)

        cov.set_xlim(0, max_gos)
        La = cov.barh((rn+2), row.mtp,
                     height=1.9, color='darkgreen',linewidth=0, alpha=0.9, label='True Positives (TP) (eggNOG)')
        Lb = cov.barh((rn+2), row.mfp, left=row.mtp,
                     height=1.9, color='indianred',linewidth=0, alpha=0.9, label='False Positives (FP) (eggNOG)')
        Lc = cov.barh((rn+2), row.munk, left=(row.mtp+row.mfp),
                     height=1.9, color='darkgrey',linewidth=0, alpha=0.9, label='Not curated terms (eggNOG)')

        Ld = cov.barh((rn-0.1), row.itp,
                     height=1.9, color='darkgreen',linewidth=0, alpha=0.4, label='(interPro)')
        Le = cov.barh((rn-0.1), row.ifp, left=row.itp,
                     height=1.9, color='indianred',linewidth=0, alpha=0.4, label='(interPro)')
        Lf = cov.barh((rn-0.1), row.iunk, left=(row.itp+row.ifp),
                     height=1.9, color='darkgrey',linewidth=0, alpha=0.4, label='(interPro)')


        iratio = row.itp / (row.itp + row.ifp + row.iunk)
        oratio = row.mtp / (row.mtp + row.mfp + row.munk)
        cov.text(max_gos, rn+2.8, '%0.2f ' %oratio, horizontalalignment='right', color = "blue")
        cov.text(max_gos, rn+0.6, '%0.2f ' %iratio, horizontalalignment='right')

        if taxa == "9606":
            cov.text(max_gos, rn+4.5, 'TP ratio\n(higher is better)',
                     horizontalalignment='right', fontweight="bold", fontsize=8)


        # ROW THIRD COLUMN
        # =========================
        cov = axes[2]
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()
        cov.set_yticks([])
        cov.set_ylim(ylim)

        cov.set_xlabel("% of proteome annotated", {'fontsize': 10})

        psize = PROTEOME_SIZE[taxa]/100
        cov.set_xlim(0, 100)

        cov.barh((rn-0.1), 100, left=0, height=4, color='grey',linewidth=0, alpha=0.05)

        Lg = cov.barh((rn+2), row.mp_tponly/psize,
                     height=1.9, color='darkblue',linewidth=0, alpha=0.9,
                     label='Proteins with TP terms only (eggNOG)')

        Lh = cov.barh((rn+2), row.mp_tpplus/psize, left=row.mp_tponly/psize,
                     height=1.9, color='purple',linewidth=0, alpha=0.9,
                     label='Proteins with TP and other terms (eggNOG)')

        Li = cov.barh((rn+2), row.mp_notp/psize, left=(row.mp_tponly + row.mp_tpplus)/psize,
                     height=1.9, color='orange',linewidth=0, alpha=0.9,
                    label='Proteins lacking TP terms (eggNOG)')

        Lj = cov.barh((rn-0.1), row.ip_tponly/psize,
                     height=1.9, color='darkblue',linewidth=0, alpha=0.4, label='(interPro)')

        Lk = cov.barh((rn-0.1), row.ip_tpplus/psize, left=row.ip_tponly/psize,
                     height=1.9, color='purple',linewidth=0, alpha=0.4, label='(interPro)')

        Lm = cov.barh((rn-0.1), row.ip_notp/psize, left=(row.ip_tponly + row.ip_tpplus)/psize,
                     height=1.9, color='orange',linewidth=0, alpha=0.4, label='(interPro)')


        iratio = row.ip_notp / (row.ip_tponly + row.ip_tpplus + row.ip_notp)
        mratio = row.mp_notp / (row.mp_tponly + row.mp_tpplus + row.mp_notp)
        cov.text(100, rn+2.8, '%0.2f ' %mratio, horizontalalignment='right', color = "blue")
        cov.text(100, rn+0.6, '%0.2f ' %iratio, horizontalalignment='right')

        if taxa == "9606":
            cov.text(100, rn+4.5, 'no-TP ratio\n(lower is better)',
                     horizontalalignment='right', fontweight="bold", fontsize=8)
        rn += 5.0

        if taxa == "9606":
            drate.legend(handles=[La, Lb, Lc, Ld, Le, Lf], loc=2, ncol=2, bbox_to_anchor=(0.1, 1.25))
            cov.legend(handles=[Lg, Lh, Li, Lj, Lk, Lm], loc=2, ncol=2, bbox_to_anchor=(-0.9, 1.25))

    plt.subplots_adjust(left=None, bottom=0, right=None, top=1.02, wspace=None, hspace=None)

    plt.savefig("plots/emapper_vs_interpro.%s.pdf" %emapper_tag, dpi=300, facecolor='w', edgecolor='w',
            orientation='portrait', bbox_inches = 'tight')

    plt.savefig("plots/emapper_vs_interpro.%s.png" %emapper_tag, dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', bbox_inches = 'tight')
    if not nodisplay:
        plt.show()

def summary_table(A, B, benchmark, cutoff):
    ''' Compare two methods and show differences in FP, TP and coverage rates '''

    tpfp_ratio_diff = []
    tp_ratio_diff = []
    terms_per_prot_diff = []
    cov_diff = []
    notp_prots_diff = []
    tponly_prots_diff = []

    for taxa in reversed(TARGET_TAXA):
        row = benchmark[taxa].loc[cutoff]

        tpfp_ratio_diff.append((row["%stp_ratio" %A] - row["%stp_ratio" %B])*100)

        tp_ratio_A = row["%stp"%A] / (row["%stp"%A] + row["%sfp"%A] + row["%sunk"%A])
        tp_ratio_B = row["%stp"%B] / (row["%stp"%B] + row["%sfp"%B] + row["%sunk"%B])
        tp_ratio_diff.append((tp_ratio_A - tp_ratio_B)*100)

        terms_A = (row["%stp"%A] + row["%sfp"%A] + row["%sunk"%A])
        terms_B = (row["%stp"%B] + row["%sfp"%B] + row["%sunk"%B])

        if terms_B > terms_A:
            terms_per_prot_diff.append((-terms_A / terms_B)*100)
        else:
            terms_per_prot_diff.append((terms_B / terms_A)*100)

        psize = PROTEOME_SIZE[taxa]/100
        cov_A = (row["%sp_tponly"%A] + row["%sp_tpplus"%A] + row["%sp_notp"%A]) /psize
        cov_B = (row["%sp_tponly"%B] + row["%sp_tpplus"%B] + row["%sp_notp"%B]) /psize
        cov_diff.append(cov_A - cov_B)

        notp_prots_A = row["%sp_notp"%A] / (row["%sp_tponly"%A] + row["%sp_tpplus"%A] + row["%sp_notp"%A])
        notp_prots_B = row["%sp_notp"%B] / (row["%sp_tponly"%B] + row["%sp_tpplus"%B] + row["%sp_notp"%B])
        notp_prots_diff.append((notp_prots_A-notp_prots_B)*100)

        tponly_prots_A = row["%sp_tponly"%A] / (row["%sp_tponly"%A] + row["%sp_tpplus"%A] + row["%sp_notp"%A])
        tponly_prots_B = row["%sp_tponly"%B] / (row["%sp_tponly"%B] + row["%sp_tpplus"%B] + row["%sp_notp"%B])
        tponly_prots_diff.append((tponly_prots_A - tponly_prots_B)*100)

    matrix = [tpfp_ratio_diff, tp_ratio_diff, terms_per_prot_diff, cov_diff, notp_prots_diff, tponly_prots_diff]
    for m in matrix:
        mean = np.mean(m)
        median = np.median(m)
        std = np.std(m)
        m.insert(0, mean)
        m.insert(1, median)
        m.insert(2, std)


    mindex = map(str.strip, """TP/FP_per_protein_ratio_diff, TP_ratio_per_protein_diff,
                               total_terms_per_prot_diff,
                               proteome_coverage_diff, notp_proteins_diff,
                               tponly_proteins_diff""".split(','))
    mcolumns=['AVERAGE', 'MEDIAN', "STD"] + [TAXA_NAMES[int(t)].split()[1] for t in TARGET_TAXA]

    return pd.DataFrame(matrix, index= mindex, columns=mcolumns)


def print_summary_table(benchmark, blast_cutoff='1e-40', target_tables=None, average_only=True):
    OUT = StringIO()
    if not target_tables:
        target_tables = set(["blast", "interpro"])

    if 'blast' in target_tables:
        print >>OUT, '\neggNOG vs BLAST (% diff)\n==============================='
        blast = summary_table("o", "h", benchmark, blast_cutoff)
        print >>OUT, blast.AVERAGE
        if not average_only:
            print >>OUT, '\neggNOG vs BLAST (% diff)\n==============================='
            print >>OUT, blast

    if 'interpro' in target_tables:
        print >>OUT, '\neggNOG vs interpro (% diff)\n============================'
        interpro = summary_table("m", "i", benchmark, blast_cutoff)
        print >>OUT, interpro.AVERAGE
        if not average_only:
            print '>>OUT, \neggNOG vs interpro (% diff)\n============================'
            print >>OUT, interpro
    return OUT

def plot_cafa2_benchmark():
    cco = [
    [0.666, 0.01, "GORBI"],
    [0.666, 0.02, "ProFun"],
    [0.612, 0.11, "PFPDB"],
    [0.489, 0.01, "APRICOT"],
    [0.471, 0.31, "eggNOG-mapper"],
    [0.467, 0.98, "EVEX"],
    [0.462, 0.98, "Tiran Lab"],
    [0.462, 0.11, "Gough Lab"],
    [0.455, 0.98, "MS-kNN"],
    [0.454, 0.89, "Orengo-FunFHMMer-mda"],
    [0.448, 0.94, "PULP"],
    [0.462, 1.00, "Naive"],
    [0.349, 0.99, "BLAST"],
    ]

    mfo = [
    [0.645, 0.25, "eggNOG-mapper"],
    [0.623, 0.05, "CBRG"],
    [0.624, 0.33, "Paccanaro"],
    [0.606, 0.22, "MS-kNN"],
    [0.606, 0.04, "SIFTER 2.4"],
    [0.605, 0.95, "Tiran Lab"],
    [0.596, 0.99, "EVEX"],
    [0.588, 0.92, "Go2Proto"],
    [0.583, 0.56, "Orengo-FunFHMMer"],
    [0.586, 0.31, "PFPDB"],
    [0.577, 0.30, "SIAM"],
    [0.348, 1.00, "Naive"],
    [0.476, 0.99, "BLAST"]
    ]

    bpo = [
    [0.461, 0.09, "GORBI"],
    [0.386, 0.27, "SIAM"],
    [0.377, 0.69, "Blast2Go"],
    [0.375, 0.96, "Tian Lab"],
    [0.374, 0.34, "eggNOG-mapper"],
    [0.372, 1.00, "Paccanaro Lab"],
    [0.366, 0.98, "MS-kNN"],
    [0.365, 0.72, "SANS"],
    [0.364, 0.93, "Orengo-FunFHMMer-mda"],
    [0.358, 0.85, "EVEX"],
    [0.352, 0.12, "ProFun"],
    [0.286, 1.00, "Naive"],
    [0.252, 0.99, "BLAST"]
    ]

    width = 0.80      # the width of the bars
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)
    fig.set_size_inches(16.9, 6)

    for target, target_name, ax in [(mfo, "Molecular Function", ax1),
                                    (bpo, "Biological Process", ax2),
                                    (cco, "Cellular Component", ax3)]:
        fmax = [v[0] for v in target]
        ind = np.arange(len(target))
        labels = [v[2] for v in target]
        ind = ind + 0.4
        colors = []
        for lab in labels:
            if lab == 'eggNOG-mapper':
                colors.append('darkgreen')
            elif lab == 'Naive':
                colors.append('darkred')
            elif lab == 'BLAST':
                colors.append('darkblue')
            else:
                colors.append('#666666')

        rects = ax.bar(ind, fmax, width, color=colors, edgecolor='none')

        # add some text for labels, title and axes ticks
        ax.set_ylabel('F-max')
        ax.set_title(target_name)
        ax.set_xticks(ind + width / 2)
        ax.set_xticklabels([v[2] for v in target])
        ax.set_ylim(0, 0.7)
        plt.sca(ax)
        plt.xticks(ind+(width/2.0)+0.1, labels, rotation='vertical')

        for i, rect in enumerate(rects):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2. - 0.2, 0.2,
                    'C=%0.2f' % target[i][1], color="white",
                    ha='left', rotation=90)

    plt.savefig("plots/emapper_cafa2.pdf", facecolor='w', edgecolor='w',
                bbox_inches = 'tight')
    plt.show()
    

def translate_tag(tag, self):
    method, orthologs, tax_scope, go_evidence = tag.split('.')
    tax_scope = 'none' if tax_scope == 'NOG' else 'auto'
    html = "<h3>%s</h3>" %tag
    html +=  "<b>emapper_method</b> = %s <br>" %method
    html += "<b>source orthologs</b> = %s <br>" %orthologs
    html += "<b>taxonomic restriction</b> = %s <br>" %tax_scope
    html += "<b>GO evidence code</b> = %s <br>" %go_evidence
    html += "<b>self proteome excluded</b> = %s <br>" %self
    return html

def translate_tag_md(tag, self):
    method, orthologs, tax_scope, go_evidence = tag.split('.')
    tax_scope = 'none' if tax_scope == 'NOG' else 'auto'
    html = "### %s\n```\n" %tag
    html += "emapper_method         = %s\n" %method
    html += "source orthologs       = %s\n" %orthologs
    html += "taxonomic restriction  = %s\n" %tax_scope
    html += "GO evidence code       = %s\n" %go_evidence
    html += "self proteome excluded = %s\n```" %self
    return html



def get_emapper_blast_summary(bench, refresh_plots=True):
    # Generate plots for all emapper runs without taxonomic restriction (for BLAST comparison purposes)
    blast_emapper_tags = sorted([k for k in bench.keys() if '.NOG.' in k])

    # Browse plots in tabs
    tabs = []
    tab_names = []
    for tag in blast_emapper_tags:
        if refresh_plots:
            plot_blast_benchmark(bench[tag], TARGET_TAXA, EVALUE_CUTOFFS, "E-value Cutoff", tag, nodisplay=True)
        stats = print_summary_table(bench[tag], '1e-40', ['blast'], average_only=True)
        html = translate_tag(tag, "yes")
        html += "<img style='width:95%%;' src='plots/emapper_vs_blast.%s.png'>" %(tag)
        html += '<br><pre>%s</pre>' %stats.getvalue()
        tabs.append(widgets.HTML(value=html))
        tab_names.append(tag)

    # Crate tab widget
    blast_tabs = widgets.Tab(children=tabs)
    for i, name in enumerate(tab_names):
        blast_tabs.set_title(i, "% 25s" %name.replace('.NOG.', '.'))
    return blast_tabs

def get_emapper_blast_summary_pdf(bench, refresh_plots=True):
    # Generate plots for all emapper runs without taxonomic restriction (for BLAST comparison purposes)
    blast_emapper_tags = sorted([k for k in bench.keys() if '.NOG.' in k])
    for tag in blast_emapper_tags:
        if refresh_plots:
            plot_blast_benchmark(bench[tag], TARGET_TAXA, EVALUE_CUTOFFS, "E-value Cutoff", tag, nodisplay=True)
        stats = print_summary_table(bench[tag], '1e-40', ['blast'], average_only=True)
        html = translate_tag_md(tag, "yes")
        html += '\n```\n %s``` ' %stats.getvalue()
        display(Markdown(html))
        display(Latex(pagebreak))
        img = '\n![%s](plots/emapper_vs_blast.%s.png)' %(tag, tag)
        display(Markdown(img))
        display(Latex(pagebreak))

def get_emapper_interpro_summary_pdf(bench, refresh_plots=True):
    interpro_emapper_tags = sorted([k for k in bench.keys() if '.auto.' in k])

    for tag in interpro_emapper_tags:
        if refresh_plots:
            plot_interpro_benchmark(bench[tag], TARGET_TAXA, tag, nodisplay=True)
        stats = print_summary_table(bench[tag], '1e-40', ['interpro'], average_only=True)
        html = translate_tag_md(tag, "yes")
        html += '\n```\n %s``` ' %stats.getvalue()
        display(Markdown(html))
        display(Latex(pagebreak))
        img = '\n![%s](plots/emapper_vs_interpro.%s.png)' %(tag, tag)
        display(Markdown(img))
        display(Latex(pagebreak))


def get_emapper_interpro_summary(bench, refresh_plots=True):
    # Generate plots for all emapper runs with default taxonomic adjustment(for interProScan comparison purposes)
    interpro_emapper_tags = sorted([k for k in bench.keys() if '.auto.' in k])

    # Browse plots in tabs
    tabs = []
    tab_names = []
    for tag in interpro_emapper_tags:
        if refresh_plots:
            plot_interpro_benchmark(bench[tag], TARGET_TAXA, tag, nodisplay=True)
        stats = print_summary_table(bench[tag], '1e-40', ['interpro'], average_only=True)
        html = translate_tag(tag, "yes")
        html += "<img style='width:95%%;' src='plots/emapper_vs_interpro.%s.png'>" %(tag)
        html += '<br><pre>%s</pre>' %stats.getvalue()
        tabs.append(widgets.HTML(value=html))
        tab_names.append(tag)

    # Crate tab widget
    interpro_tabs = widgets.Tab(children=tabs)
    for i, name in enumerate(tab_names):
        interpro_tabs.set_title(i, name)

    return interpro_tabs

def plot_general_benchmark(benchmark, target_taxa, emapper_tag):
    fig_height = len(target_taxa) + 1
    f, axes = plt.subplots(1, 2, figsize=(16.9, fig_height))
    #max_gos = max([(benchmark[t].htp + benchmark[t].hfp + benchmark[t].hunk)[0] for t in target_taxa])
    max_gos = 0
    for taxa in target_taxa:
        row = benchmark[taxa].loc["1e-40"]
        max_gos = max(max_gos, row.htp + row.hfp + row.hunk)

    rn = 0
    labels = []
    for taxa in reversed(target_taxa):
        labels.append(' '.join(TAXA_NAMES.get(int(taxa), taxa).split()[:2]))

        row = benchmark[taxa].loc["1e-40"]
        yticks = [3]
        for t in range(len(target_taxa)-1):
            yticks.append(yticks[-1] + 1 + 6)
        ylim = (0, len(target_taxa)*7)


        # ROW IN SECOND COLUMN
        # =========================

        cov = axes[0]
        cov.set_xlabel("Average GO terms per protein", {'fontsize': 10})
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()
        cov.set_yticks([])
        cov.set_ylim(ylim)

        cov.set_xlim(0, max_gos)
        La = cov.barh((rn+4), row.mtp,
                     height=1.9, color='darkgreen',linewidth=0, alpha=0.9, label='True Positives (TP)')
        Lb = cov.barh((rn+4), row.mfp, left=row.mtp,
                     height=1.9, color='indianred',linewidth=0, alpha=0.9, label='False Positives (FP)')
        Lc = cov.barh((rn+4), row.munk, left=(row.mtp+row.mfp),
                     height=1.9, color='darkgrey',linewidth=0, alpha=0.9, label='Not curated terms')

        Ld = cov.barh((rn+1.9), row.itp,
                     height=1.9, color='darkgreen',linewidth=0, alpha=0.4)
        Le = cov.barh((rn+1.9), row.ifp, left=row.itp,
                     height=1.9, color='indianred',linewidth=0, alpha=0.4)
        Lf = cov.barh((rn+1.9), row.iunk, left=(row.itp+row.ifp),
                     height=1.9, color='darkgrey',linewidth=0, alpha=0.4)

        Lg = cov.barh((rn-0.2), row.htp,
                      height=1.9, color='darkgreen',linewidth=0, alpha=0.4)
        Lh = cov.barh((rn-0.2), row.hfp, left=row.htp,
                      height=1.9, color='indianred',linewidth=0, alpha=0.4)
        Li = cov.barh((rn-0.2), row.hunk, left=(row.htp+row.hfp),
                      height=1.9, color='darkgrey',linewidth=0, alpha=0.4)


        iratio = row.itp / (row.itp + row.ifp + row.iunk)
        oratio = row.mtp / (row.mtp + row.mfp + row.munk)
        hratio = row.htp / (row.htp + row.hfp + row.hunk)

        cov.text(max_gos, rn+4.8, 'eMapper        %0.2f ' %oratio, horizontalalignment='right', color = "blue")
        cov.text(max_gos, rn+2.6, 'interPro       %0.2f ' %iratio, horizontalalignment='right')
        cov.text(max_gos, rn+0.5, 'BLAST 10E-40   %0.2f ' %hratio, horizontalalignment='right')

        
        if taxa == target_taxa[0]:
            cov.text(max_gos, rn+6.5, 'TP ratio\n(higher is better)',
                     horizontalalignment='right', fontweight="bold", fontsize=8)
            #drate.get_yaxis().set_visible(False)
            cov.spines['right'].set_visible(False)
            cov.spines['top'].set_visible(False)
            cov.spines['left'].set_visible(False)
            cov.get_xaxis().tick_bottom()
            cov.get_yaxis().tick_left()
            cov.set_xlabel("Average GO terms per protein", {'fontsize': 10})
            cov.set_yticks(yticks) #[2, 7, 12, 17, 22])
            cov.set_yticklabels(labels, fontsize=12, style='italic')

            cov.set_ylim(ylim)
            #cov.yaxis.tick_right()


        # ROW THIRD COLUMN
        # =========================
        cov = axes[1]
        cov.spines['right'].set_visible(False)
        cov.spines['top'].set_visible(False)
        cov.spines['left'].set_visible(False)
        cov.get_xaxis().tick_bottom()
        cov.set_yticks([])
        cov.set_ylim(ylim)

        cov.set_xlabel("% of proteome annotated", {'fontsize': 10})

        psize = PROTEOME_SIZE[taxa]/100
        cov.set_xlim(0, 100)

        #cov.barh((rn-0.1), 100, left=0, height=4, color='grey',linewidth=0, alpha=0.05)

        Lg = cov.barh((rn+4), row.mp_tponly/psize,
                     height=1.9, color='darkblue',linewidth=0, alpha=0.9,
                     label='Proteins with TP terms only')

        Lh = cov.barh((rn+4), row.mp_tpplus/psize, left=row.mp_tponly/psize,
                     height=1.9, color='purple',linewidth=0, alpha=0.9,
                     label='Proteins with TP and other terms')

        Li = cov.barh((rn+4), row.mp_notp/psize, left=(row.mp_tponly + row.mp_tpplus)/psize,
                     height=1.9, color='orange',linewidth=0, alpha=0.9,
                    label='Proteins lacking TP terms')

        Lj = cov.barh((rn+1.9), row.ip_tponly/psize,
                     height=1.9, color='darkblue',linewidth=0, alpha=0.4)

        Lk = cov.barh((rn+1.9), row.ip_tpplus/psize, left=row.ip_tponly/psize,
                     height=1.9, color='purple',linewidth=0, alpha=0.4)

        Lm = cov.barh((rn+1.9), row.ip_notp/psize, left=(row.ip_tponly + row.ip_tpplus)/psize,
                     height=1.9, color='orange',linewidth=0, alpha=0.4)


        Ln = cov.barh((rn-0.2), row.hp_tponly/psize,
                      height=1.9, color='darkblue',linewidth=0, alpha=0.4)

        Lo = cov.barh((rn-0.2), row.hp_tpplus/psize, left=row.hp_tponly/psize,
                      height=1.9, color='purple',linewidth=0, alpha=0.4)

        Lp = cov.barh((rn-0.2), row.hp_notp/psize, left=(row.hp_tponly + row.hp_tpplus)/psize,
                      height=1.9, color='orange',linewidth=0, alpha=0.4)


        iratio = row.ip_notp / (row.ip_tponly + row.ip_tpplus + row.ip_notp)
        mratio = row.mp_notp / (row.mp_tponly + row.mp_tpplus + row.mp_notp)
        hratio = row.hp_notp / (row.hp_tponly + row.hp_tpplus + row.hp_notp)
        cov.text(100, rn+4.8, 'eMapper         %0.2f ' %mratio, horizontalalignment='right', color = "blue")
        cov.text(100, rn+2.6, 'interPro        %0.2f ' %iratio, horizontalalignment='right')
        cov.text(100, rn+0.5, 'BLAST 10E-40    %0.2f ' %hratio, horizontalalignment='right')

        if taxa == target_taxa[0]:
            cov.text(100, rn+6.5, 'no-TP ratio\n(lower is better)',
                     horizontalalignment='right', fontweight="bold", fontsize=8)
        rn += 7.0

        if taxa == target_taxa[0]:
            cov.legend(handles=[La, Lb, Lc, Lg, Lh, Li], loc=2, ncol=2, bbox_to_anchor=(-0.9, 1.50))

    plt.subplots_adjust(left=None, bottom=0, right=None, top=1.02, wspace=None, hspace=None)
    plt.savefig("plots/emapper_non_model_species.%s.pdf" %(emapper_tag), facecolor='w', edgecolor='w',
            orientation='portrait', bbox_inches = 'tight')

    plt.show()

