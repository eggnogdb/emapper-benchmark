import sys
from eggnogmapper import annota
annota.connect()
GO_EXCLUDED =  set(["ND", "IEA"])
GO_EXPERIMENTAL = set(["EXP","IDA","IPI","IMP","IGI","IEP"])

for species in map(int, sys.argv[1:]):
    print 'Processing', species
    proteome_file = "data/proteomes/%s.fa" %species

    EXP_OUT = open('groundtruth/%s.groundtruth.experimental' %species, 'w')
    NOIEA_OUT = open('groundtruth/%s.groundtruth.non-electronic' %species, 'w')

    for line in open(proteome_file):
        if line.strip().startswith('>'):
            query = line.split()[0][1:]
            _, exp, _ = annota.get_member_annotations([query], target_go_ev=GO_EXPERIMENTAL,
                                                         excluded_go_ev=GO_EXCLUDED)
            _, noIEA, _ = annota.get_member_annotations([query], target_go_ev=None,
                                                         excluded_go_ev=GO_EXCLUDED)

            print >>EXP_OUT, '\t'.join([query, ','.join(sorted(exp))])
            print >>NOIEA_OUT, '\t'.join([query, ','.join(sorted(noIEA))])
    EXP_OUT.close()
    NOIEA_OUT.close()
