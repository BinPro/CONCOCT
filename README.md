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
###Fundamental dependencies###
```
python v2.7.*
gcc
gsl
```

These items are prerequisities for the installation of concoct as described below. The installation procedure varies on different systems, and described in this README is only how to proceed with a linux (ubuntu) distribution.

The first item, ```python v2.7.*```, should be installed on a modern Ubuntu distribution. A c-compiler, e.g. ```gcc```, is needed to compile the c parts of concoct that uses the GNU Scientific Library ```gsl```. For linux (ubuntu) this is installed through:
```
apt-get install build-essential libgsl0-dev
```
###Python packages###
```
cython>=0.19.2
numpy>=1.7.1
scipy>=0.12.0
pandas>=0.11.0
biopython>=1.62b
scikit-learn>=0.13.1
```
These are the python packages that need to be installed in order to run concoct. If you follow the installation instructions below, these will be installed automatically, but are listed here for transparency. 

###Optional dependencies###

* To create the input table (containing average coverage per sample and contig)
    * [BEDTools](https://github.com/arq5x/bedtools2/releases) version >= 2.15.0 (only genomeCoverageBed)
    * [Picard](https://launchpad.net/ubuntu/+source/picard-tools/) tools version >= 1.77
    * [samtools](http://samtools.sourceforge.net/) version >= 0.1.18
    * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) version >= 2.1.0
    * [GNU parallel](http://www.gnu.org/software/parallel/) version >= 20130422

* For validation of clustering using single-copy core genes
    * [PROKKA](http://www.vicbioinformatics.com/software.prokka.shtml)
    * Python packages: ```bcbio-gff>=0.4```
    * R packages: ```gplots, reshape, ggplot2, ellipse, getopt``` and ```grid```

##Installation#
Here we describe two recommended ways of getting concoct to run on your computer/server. The first option, using Anaconda, should work for any *nix (e.g. Mac OS X or Linux) system even where you do not have 'sudo' rights (e.g. on a common computer cluster). The second option is suitable for a linux computer where you have root privileges and you prefer to use a virtual machine where all dependencies to run concoct are included.

###Using Anaconda###
This instruction shows how to install all dependencies (except the 'Fundamental dependencies' and the 'Optional dependencies' listed above) using an Anaconda environment. Anaconda is a tool to isolate your python installation, which allows you to have multiple parallel installations using different versions of different packages, and gives you a very convenient and fast way to install the most common scientific python packages. Anaconda is free but not open source, you can download Anaconda [here](https://store.continuum.io/cshop/anaconda/). Installation instructions can be found [here](http://docs.continuum.io/anaconda/install.html).

After installing Anaconda, create a new environment that will contain the concoct installation:
```
conda create -n concoct_env python=2.7.6
```
After choosing to proceed, run the suggested command:
```
source activate concoct_env
```
then install the concoct dependencies into this environment:
```
conda install cython numpy scipy biopython pandas pip scikit-learn
```
Finally, download the CONCOCT distribution from https://github.com/BinPro/CONCOCT/releases (stable) and extract the files, or clone the repository with github (potentially unstable). Resolve all dependencies, see above and then execute within the CONCOCT directory:
```
python setup.py install
```

###Using Docker###
If you have root access where you want to install concoct and storage for roughly 1.2G "virtual machine" then Docker provides a very nice way to get a Docker image with concoct and its dependencies installed. This way the only thing you install on your host system is Docker, the rest is contained in an Docker image. This allows you to install and run programs in that image without it affecting your host system. You should get to know Docker here: https://www.docker.io/the_whole_story/
You need to get Docker installed (see https://www.docker.io/gettingstarted/ and specially if you have Ubuntu http://docs.docker.io/en/latest/installation/ubuntulinux/). When Docker is installed you need to download and log into the concoct image which can be done in one command. We also want to map a folder from the host (/home/user/MyData) to a folder in the image (/opt/MyData). To get all this working we execute one command:
```
sudo docker run -v /home/user/MyData:/opt/MyData -i -t binnisb/concoct_0.2.1 bash
```
This downloads the image (about 1.2G) and logs you into a bash shell. To test concoct you can then do:
```
$ cd /opt/CONCOCT-0.2.1
$ nosetests
```
Which should execute all tests without errors. Then to run concoct on your data (stored in /home/user/MyData on host) you can do:
```
$ cd /opt/MyData
$ concoct --coverage_file coverage.csv --composition_file composition.fa -b output_folder/
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

