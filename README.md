Scripts and analysis of [eggNOG-mapper](https://github.com/jhcepas/eggnog-mapper) performance compared with BLAST and InterProScan.
Supplementary Material for the paper: "Fast genome-wide functional annotation through orthology assignment by eggNOG-mapper" http://biorxiv.org/content/early/2016/09/22/076331

# Analysis

Additional figures and details are available as [jupyter notebook](./benchmark_analysis.ipynb) in this reposotory. 


# Reproducibility

Benchmark requires running BLAST, InterProScan and eggNOG-mapper independently
and, in the case of eggNOG-mapper, using multiple combinations of parameters.
All raw data (~55GB uncompressed) can be downloaded from:
http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/emapper_benchmark_data.tar.gz

Data includes 4 self-explanatory directories: 

```bash
  blast/
  interpro/
  emapper_hmm/
  emapper_diamond/
```

Each directory contains a `run.sh` script with the command lines used to generate all data. 

To re-run the whole benchmark analysis, decompress all 4 directories under
within this repository, and execute `benchmark.py`, which will generate a
digested `all_benchmark_tables.pkl` file to be used within the jupyter notebook.





