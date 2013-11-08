#CONCOCT [![Build Status](https://travis-ci.org/BinPro/CONCOCT.png?branch=master)](https://travis-ci.org/BinPro/CONCOCT)#

A program for unsupervised binning of metagenomic contigs by using nucleotide composition, 
coverage data in multiple samples and linkage data from paired end reads.

Warning! This software is to be considered unstable and under heavy development. Functionality can not yet be guaranteed and the user interface may change significantly. 
If you still want to use this software, please stay up to date with the list of known issues:
https://github.com/BinPro/CONCOCT/issues

##Install##
###Short version###
Installs the package concoct in default python path, and adds script concoct to bin. You can use sudo if needed.

Resolve all dependencies, see below and then clone the repository and execute:
```
cd CONCOCT
python setup.py install
```

###Detailed version###
Installs the package concoct in default python path, and adds script concoct to bin. You can use sudo if needed.

The simplest way to get the dependencies (given Ubuntu / Debian, similar for other distros):
```
sudo apt-get install git python-setuptools python-biopython python-nose \
                     python-numpy python-pandas python-scikits-learn python-scipy
git clone https://github.com/BinPro/CONCOCT.git
cd CONCOCT
python setup.py install
```
And if you want MPI execution:
```
sudo apt-get install python-pip openmpi1.6-bin libopenmpi1.6 libopenmpi1.6-dev
sudo pip install mpi4py
```

##Execute concoct##
The script concoct takes two input files. The first file, the coverage
file, contains a table where each row correspond to a contig, and each
column correspond to a sample. The values are the average coverage for
this contig in that sample. All values are separated with tabs. The second file contains sequences in fasta format. It is named the 
composition file since it is used to calculate the kmer composition,
or the genomic signature, of each contig.

Here is a list of all parameters available for the concoct script. For a complete and up to date 
explanation of each parameter and option, the recommended way is to run

```
usage: concoct [-h] [--coverage_file COVERAGE_FILE]
               [--composition_file COMPOSITION_FILE] [-c CLUSTERS]
               [-k KMER_LENGTH] [-l LENGTH_THRESHOLD] [-r READ_LENGTH] [-s]
               [--coverage_percentage_pca COVERAGE_PERCENTAGE_PCA]
               [--composition_percentage_pca COMPOSITION_PERCENTAGE_PCA]
               [--total_percentage_pca TOTAL_PERCENTAGE_PCA] [-e EXECUTIONS]
               [-i ITERATIONS] [-b BASENAME] [-p] [-m MAX_N_PROCESSORS]
               [-f FORCE_SEED] [--covariance_type {full,diag}]
               [--no_cov_normalization] [--no_total_coverage]
```

For a complete and up to date explanation of each parameter and option, the recommended way is to run


```
concoct --help
```

##Dependencies##

concoct requires python version 2.7 and the following packages:
```
argparse==1.2.1
biopython==1.62b
nose==1.3.0
numpy==1.7.1
pandas==0.11.0
scikit-learn==0.13.1
scipy==0.12.0
```
If mpi will be used for parallelization, also add the python package <pre>mpi4py==1.3.1</pre> and linux (ubuntu) repositories:
```
openmpi1.6-bin 
libopenmpi1.6-dev
```

