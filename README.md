CONCOCT
======
A program for unsupervised binning of metagenomic contigs by using nucleotide composition, 
coverage data in multiple samples and linkage data from paired end reads.

Warning! This software is to be considered unstable and under heavy development. Functionality can not yet be guaranteed and the user interface may change significantly. 
If you still want to use this software, please stay up to date with the list of known issues:
https://github.com/BinPro/CONCOCT/issues

Install
-------
Clone the repository and execute
```
cd CONCOCT
python setup.py install
```
Installs the package concoct in default python path, and adds script CONCOCT to bin

Execute concoct
-------
```
positional arguments:
  coverage_file   specify the coverage file
  composition_file  specify the composition file
optional arguments:
  -h, --help                        show help message and exit
  -c --clusters                     specify range of clusters to try out on format
                                    first,last,step. default 20,100,2.
  -n --coverage_file_column_names   specify the first and last column names for continuous
                                    coverage range of read counts as first,last
  -k --kmer_length                  specify kmer length, defaults to tetramer
  -l --limit_kmer_count             specify the kmer count for threshold in running PCA on
                                    composition contigs, default 1000
  -r --read_length                  specify read length for coverage, default 100
  -s --split_pca                    specify this flag to first do PCA for the composition and
                                    using that component number that explaines 90 percent of
                                    variance for the coverage as well. Default join
                                    composition and coverage before PCA. ***NOT IMPLEMENTED***
  -e --executions                   How often to initialize each cluster count. default 5
                                    times
  -i --iterations                   Maximum number of iterations if convergance not achieved
  -b --basename                     Specify the basename for files or directory where output
                                    will be placed. Path to existing directory or basename
                                    with a trailing '/' will be interpreted as a directory.
                                    If not provided, current directory will be used.
  -p --pipe                         Add this tag if the main result file should be
                                    printed to stdout. Useful for pipeline use.
  -m --max_n_processors             Specify the maximum number of processors CONCOCT is
                                    allowed to use, if absent, all present processors will
                                    be used.
  -f --force_seed                   Specify an integer to use as seed for clustering to
                                    make the process deterministic. Default is 0 which
                                    makes the process random.
```

Dependencies
-----------
CONCOCT requires python version 2.7 and the following packages:
```
argparse==1.2.1
biopython==1.62b
nose==1.3.0
numpy==1.7.1
pandas==0.11.0
scikit-learn==0.13.1
scipy==0.12.0
```
