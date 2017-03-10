
This supplementary data accompanies the results of the CAFA2 experiment, as described in "An expanded evaluation of protein function prediction methods shows an improvement in accuracy" (Jiang et al, submitted).

The aim is to provide (1) transparency (2) reproducibility and (3) archiving of the data used for the
CAFA2 experiment. The data are provided under the CC-0 licence.

The following describes the data in the repository. 

Many of the file contents are described by their names, using the following abbreviations which are elaborated upon in the paper:

MFO: Molecular Function Ontology
BPO: Biological Process Ontology
CCO: Cellular Component Ontology
NK: no-knowledge targets
LK: limited-knowledge targets
full: full mode (all) targets
partial: partial mode (>5,000 targets)
auc: area under the curve, term centric analysis
fmax: analysis scored by Fmax.
wfmax: analysis scored by weighted Fmax.
nsmin: analysis scored by normalized Smin
smin: analysis scored by Smin

Category code:
all: all benchmark
eukarya: all eukaryotic organisms
prokarya: all prokaryotic organisms
easy: all easy targets
hard: all difficult targets
ARATH: Arabidopsis thaliana
DANRE: Danio rerio
DICDI: Dictyostelium discoideum
DROME: Drosophila melanogaster
ECOLI: Escherichia coli K12
HUMAN: Homo sapiens
MOUSE: Mus musculus
PSEAE: Pseudomonas aeruginosa
RAT: Rattus norvegicus
SCHPO: Schizosaccharomyces pombe
YEAST: Saccharomyces cerevisiae


Therefore, for example, cco_RAT_LK_full_all_smin_sheet.csv, is a table describing the results for the evaluation of the targets in Rat, cellular component ontology, limited-knowledge, full mode, Smin used for scoring.  


Directory contents:
.
|-- data
|   |-- benchmark lists of proteins for each benchmark category
|   |-- CAFA2-targets lists of CAFA2 target proteins
|   |-- GO-t0 parsed annotations from GO Consortium used for CAFA at submission time
|   |-- GO-t1 parsed annotations from GO Consortium used for CAFA at scoring time
|   |-- HPO-t0 Human Phenotype Ontology annotations used for CAFA at submission time
|   |-- HPO-t1 Human Phenotype Ontology annotations used for CAFA at scoring time
|   |-- ontology The four ontologies used in CAFA experiment
|   |-- SwissProt-t0 SwissProt experimental annotations used for CAFA at submission time
|   |-- SwissProt-t1 SwissProt experimental annotations used for CAFA at scoring time
|   |-- UniProt-GOA-t0 UniProt-GOA annotations used for CAFA at submission time
|   `-- UniProt-GOA-t1 UniProt-GOA annotations used for CAFA at scoring time
|-- plot all result plots
`-- sheet all result tables

CAFA code can be found in Github:
https://github.com/yuxjiang/CAFA2

