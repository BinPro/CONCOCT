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
TODO! Not correct flags !!!!!!!!!!!!!!!!!!!!!!!!!
- f               contains the fasta formatted contigs
- k               the kmer size
- c               number of clusters
- r               the number of runs to execute clustering
- i               the number of maximum iterations per run of clustering
- cf              file with the coverage data
- o               The directory where result files should be stored, otherwise current dir used
```

Dependencies
-----------
Developed under python 2.7 and following packages installed through pip:
```
numpy==1.7.1
pandas==0.11.0
biopython==1.61
scikit-learn
scipy==0.12.0
argparse==1.2.1
```
