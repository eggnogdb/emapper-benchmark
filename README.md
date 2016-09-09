Scripts and analysis of [eggNOG-mapper](https://github.com/jhcepas/eggnog-mapper) performance compared with BLAST and InterProScan. 

# Analysis

Additional figures and details are available as [jupyter notebook](./benchmark_analysis.ipynb) in this reposotory. 


# Reproducibility

Benchmark worflow consist of the following steps (not necessarily in the same order and with the same syntax as shown here):

### Run blast searches

```bash
for x in 9606 3702 7227 4932 511145; do    
    blastp -db eggnog4.proteins.core_periphery.fa -query $x.fa -evalue 0.001 -num_threads 20 -outfmt 6 -out $x.hits -max_target_seqs 1000000;
done;
```

### Run interpro annotations

```bash
for x in 9606 3702 7227 4932 511145; do    
    interproscan.sh -i $x.fa -f TSV --goterms --iprlookup -pa -o $x.interpro;
done;
```

### Run eggnog-Mapper annotations


```bash
for x in 9606 3702 7227 4932; do    
    python eggnog-mapper/emapper.py -i $x.fa --cpu 20 --output $x.emapper -d euk;
done

for x in 511145; do    
    python eggnog-mapper/emapper.py -i $x.fa --cpu 20 --output $x.emapper -d bact;
done

```

### Run `compute_blast_benchmark.py` to generate the benchmark files under
  `benchmark_results`, where different cutoffs are applied to blast and
  eggnog-mapper results.

```bash
for x in 9606 3702 7227 4932 511145; do    
    python compute_blast_benchmark.py $x.hits $x.emapper.refined_hits.annot $x identity;        
    python compute_blast_benchmark.py $x.hits $x.emapper.refined_hits.annot $x evalue;
done;
```


###Run `pool_[emapper|interpro]_annotations.py`

... to extract raw annotations and augment the list of GO terms for each query
  so parent terms are always included (otherwise interpro would look like
  producing way less annotations).

```bash
python pool_interpro_annotations.py results_1.0.0
python pool_interpro_annotations.py results_1.0.0
```




