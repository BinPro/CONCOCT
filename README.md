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
usage: CONCOCT [-h] [-c CLUSTERS] [-n COVERAGE_FILE_COLUMN_NAMES]
               [-k KMER_LENGTH] [-l LIMIT_KMER_COUNT] [-r READ_LENGTH] [-s]
               [--coverage_percentage_pca COVERAGE_PERCENTAGE_PCA]
               [--composition_percentage_pca COMPOSITION_PERCENTAGE_PCA]
               [--total_percentage_pca TOTAL_PERCENTAGE_PCA] [-e EXECUTIONS]
               [-i ITERATIONS] [-b BASENAME] [-p] [-m MAX_N_PROCESSORS]
               [-f FORCE_SEED]
               coverage_file composition_file

positional arguments:
  coverage_file         specify the coverage file
  composition_file      specify the composition file

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        specify range of clusters to try out on format
                        first,last,step. default 20,100,2.
  -n COVERAGE_FILE_COLUMN_NAMES, --coverage_file_column_names COVERAGE_FILE_COLUMN_NAMES
                        specify the first and last column names for continuous
                        coverage range of read counts as first,last
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        specify kmer length, defaults to tetramer
  -l LIMIT_KMER_COUNT, --limit_kmer_count LIMIT_KMER_COUNT
                        specify the kmer count for threshold in running PCA on
                        composition contigs, default 1000
  -r READ_LENGTH, --read_length READ_LENGTH
                        specify read length for coverage, default 100
  -s, --split_pca       specify this flag to perform PCA for the composition
                        and coverage data seperately
  --coverage_percentage_pca COVERAGE_PERCENTAGE_PCA
                        The percentatage of variance explained by the
                        principal components for the coverage data, only
                        considered if split pca is used
  --composition_percentage_pca COMPOSITION_PERCENTAGE_PCA
                        The percentatage of variance explained by the
                        principal components for the composiiton data, only
                        considered if split pca is used
  --total_percentage_pca TOTAL_PERCENTAGE_PCA
                        The percentage of variance explained by the principal
                        components for the combined data, only considered if
                        split pca is NOT used
  -e EXECUTIONS, --executions EXECUTIONS
                        How often to initialize each cluster count. default 5
                        times
  -i ITERATIONS, --iterations ITERATIONS
                        Maximum number of iterations if convergance not
                        achieved
  -b BASENAME, --basename BASENAME
                        Specify the basename for files or directory where
                        outputwill be placed. Path to existing directory or
                        basenamewith a trailing '/' will be interpreted as a
                        directory.If not provided, current directory will be
                        used.
  -p, --pipe            Add this tag if the main result file should beprinted
                        to stdout. Useful for pipeline use
  -m MAX_N_PROCESSORS, --max_n_processors MAX_N_PROCESSORS
                        Specify the maximum number of processors CONCOCT is
                        allowed to use, if absent, all present processors will
                        be used. Default is to assume MPI execution, allowing
                        execution over multiple nodes, with a fallback to
                        standard python multiprocessing on a single machine
                        using available cores. It is recommended to install
                        mpi and mpi4py if execution over multiple nodes is
                        required. To run with MPI use call it with mpirun -np
                        N or equivalent
  -f FORCE_SEED, --force_seed FORCE_SEED
                        Specify an integer to use as seed for clustering. You
                        can specify 0 for random seed. The default seed is 11.
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
mpi4py==1.3.1
```
If mpi will be used for parallelization, also add the python package <pre> mpi4py==1.3.1 </pre> and linux (ubuntu) repositories:
```
openmpi1.6-bin 
libopenmpi1.6-dev
```

