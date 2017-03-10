# Based on direct mapping of reads to original genes (by genomic coordinates) in
# the reference genomes --> Calculate gene simulate gene abundances in each
# sample
for x in `seq 0 3`; do
    echo $x; zcat sample-$x.read_genes.gz|awk '{if ($5=="True" && $4=="True"){print $line}}'|cut -f3|sort|datamash -g 1 count 1 > sample-$x.original_gene_abundances;
done
