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
TODO! Not correct dependencies!!!!!!!!!!!!!!
Jinja2==2.6
Pygments==1.6
argparse==1.2.1
biopython==1.61
distribute==0.6.35
docutils==0.10
ipython==0.13.2
line-profiler==1.0b3
matplotlib==1.2.1
nose==1.3.0
numpy==1.7.1
openpyxl==1.6.2
pandas==0.11.0
python-dateutil==2.1
pytz==2013b
pyzmq==13.1.0
scipy==0.12.0
six==1.3.0
stevedore==0.8
tornado==3.0.1
virtualenv==1.9.1
virtualenv-clone==0.2.4
virtualenvwrapper==3.7
wsgiref==0.1.2
```
