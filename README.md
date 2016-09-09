Scripts and analysis of [eggNOG-mapper](https://github.com/jhcepas/eggnog-mapper) performance compared with BLAST and InterProScan. 

# Analysis

Additional figures and details are available as [jupyter notebook](./benchmark_analysis.ipynb) in this reposotory. 


# Reproducibility

Benchmark requires running BLAST, InterProScan and eggNOG-mapper independently. Programs were executed as followed. 
Raw_data (~55GB uncompressed) can be downloaded from: http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/emapper_benchmark.tar.gz
The downloaded data should include a `run.sh` script with the command lines used to call each program. 

To rerun the whole benchmark analysis, decompress all the raw data in the root directory of this repository, and call `benchmark.py`, which will generate the digested `all_benchmark_tables.pkl`
file used in the jupyter notebook. 



