#CONCOCT 0.2.1 [![Build Status](https://travis-ci.org/BinPro/CONCOCT.png?branch=master)](https://travis-ci.org/BinPro/CONCOCT)#

A program for unsupervised binning of metagenomic contigs by using nucleotide composition, 
coverage data in multiple samples and linkage data from paired end reads.

Warning! This software is to be considered under development. Functionality and the user interface may still change significantly from one version to another.
If you want to use this software, please stay up to date with the list of known issues:
https://github.com/BinPro/CONCOCT/issues

##Mailing List##
Feel free to contact our mailing list concoct-support@lists.sourceforge.net for questions regarding the installation procedure and/or the usage of concoct that you feel is not sufficiently described in the online documentation. 

If you would like subscribe to concoct-support mailing list, you can do so [here](https://lists.sourceforge.net/lists/listinfo/concoct-support)

##Dependencies##

In order to install concoct, it requires python version 2.7.* and the python package installer ```pip```. It also requires a c compiler, e.g. ```gcc``` and the GNU Scientific Library ```gsl```. For linux (ubuntu) this is installed through:
```
apt-get install build-essential gsl-bin
```

Before or during the installation of concoct, several other python packages will be downloaded and installed by pip.

##Install##
Install the package concoct in default python path, and adds script concoct to bin. You can use sudo if needed.

###Using Docker and Anaconda###
Docker provides a very nice way to get a virtual machine with concoct and all dependencies installed. You need to get Docker installed (see https://www.docker.io/gettingstarted/ and specially if you have Ubuntu http://docs.docker.io/en/latest/installation/ubuntulinux/). When Docker is installed you need to download and log into the virtual machine which can be done in one command. We also want to map a folder from the host (/home/user/MyData) to a folder in the virtual machine (/opt/MyData). To get all this working we execute one command:
```
sudo docker run -v /home/user/MyData:/opt/MyData -i -t binnisb/concoct_0.2.1 bash
```
This downloads the VM (about 1.2G) and logs you into a bash shell. To test concoct you can then do:
```
$ cd /opt/CONCOCT
$ nosetests
```
Which should execute all tests without errors. Then to run concoct on your data (stored in /home/user/MyData on host) you can do:
```
$ cd /opt/MyData
$ concoct --coverage_file coverage.csv --composition_file composition.fa -b output_folder/
```

###Using Ubuntu and Anaconda###
On Ubuntu this will install all dependencies and the anaconda environment
```
apt-get update -qq
apt-get install -qq wget git build-essential libgsl0-dev gsl-bin
wget http://repo.continuum.io/miniconda/Miniconda-3.3.0-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -p /home/travis/miniconda -b
export PATH=/opt/miniconda/bin:$PATH
conda update --yes conda
conda install --yes python=2.7 atlas cython numpy scipy biopython pandas pip scikit-learn
```

###Using pip###
Download the CONCOCT distribution from https://github.com/BinPro/CONCOCT/releases (stable) and extract the files, or clone the repository with github (potentially unstable)
```
git clone https://github.com/BinPro/CONCOCT.git
```

Resolve all dependencies, see above and then execute:
```
cd CONCOCT
pip install -r requirements.txt
python setup.py install
```

###Using apt-get###
Another way to get the dependencies (given Ubuntu / Debian, similar for other distros) is through ```apt-get```. However, for some packages, only deprecated versions are available. Make sure that the requirements for these packages are fulfilled:

    biopython>=1.62b
    numpy>=1.7.1
    pandas>=0.11.0
    scikit-learn>=0.13.1
    scipy>=0.12.0

The actual commands for installing is then
```
sudo apt-get install git python-setuptools python-biopython python-nose \
                     python-numpy python-pandas python-scikits-learn python-scipy \
                     build-essential gsl-bin
git clone https://github.com/BinPro/CONCOCT.git
cd CONCOCT
python setup.py install
```

##Execute concoct##
The script concoct takes two input files. The first file, the coverage
file, contains a table where each row correspond to a contig, and each
column correspond to a sample. The values are the average coverage for
this contig in that sample. All values are separated with tabs. The second file contains sequences in fasta format. It is named the 
composition file since it is used to calculate the kmer composition,
or the genomic signature, of each contig.

Here is a list of all parameters available for the concoct script.
```
usage: concoct [-h] [--coverage_file COVERAGE_FILE]
               [--composition_file COMPOSITION_FILE] [-c CLUSTERS]
               [-k KMER_LENGTH] [-l LENGTH_THRESHOLD] [-r READ_LENGTH]
               [--total_percentage_pca TOTAL_PERCENTAGE_PCA] [-b BASENAME]
               [-s SEED] [-i ITERATIONS] [-e EPSILON] [--no_cov_normalization]
               [--no_total_coverage] [-o] [-d] [-v]
```

For a complete explanation of each parameter and option, the recommended way is to run


```
concoct --help
```

### Complete Example ###
After having installed concoct, a complete workflow can be found [here](https://github.com/BinPro/CONCOCT/blob/master/doc/complete_example.md).

