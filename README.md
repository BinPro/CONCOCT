CONCOCT
======

A program for unsupervised binning of metagenomic contigs by using nucleotide composition, 
coverage data in multiple samples and linkage data from paired end reads.

Install
-------
Clone the repository and execute
```
cd CONCOCT/src
python setup.py install
```
Installs the package concoct in default python path, and adds script concoct to bin

Execute concoct
-------
```
positional arguments:
  coverage_file   specify the coverage file
  composition_file  specify the composition file
  n                 specify the first and last column names for continuous
                    coverage range of read counts as first,last
optional arguments:
  -h, --help        show help message and exit
  -c C              specify range of clusters to try out on format
                    first,last,step. default 20,100,2.
  -k K              specify kmer length, defaults to tetramer
  -l L              specify the kmer count for threshold in running PCA on
                    composition contigs, default 1000
  -r R              specify read length for coverage, default 100
  -s                specify this flag to first do PCA for the composition and
                    using that component number that explaines 90 percent of
                    variance for the coverage as well. Default join
                    composition and coverage before PCA. ***NOT IMPLEMENTED***
  -e E              How often to initialize each cluster count. default 5
                    times
  -i I              Maximum number of iterations if convergance not achieved
  -o O              specify where output directory will be placed.
     		        If not provided, current directory will be used.
		    	    All files will be created in:
                    folder/CONCOCT_YYMMDD_HHMM_SSXXXXXX
```

Dependencies
-----------
Developed under python 2.7 and following packages installed through pip:
```
argparse==1.2.1
biopython==1.62b
nose==1.3.0
numpy==1.7.1
pandas==0.11.0
scikit-learn==0.13.1
scipy==0.12.0
```
