# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v{{version}}
# VERSION {{version}}
# 
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v{{version}} installed.
#
# To login to the docker with a shared directory from the host do:
#
# sudo docker run -v /my/host/shared/directory:/my/docker/location -i -t binnisb/concoct_{{version}} /bin/bash
#

FROM ubuntu:13.10
MAINTAINER CONCOCT developer group, concoct-support@lists.sourceforge.net

ENV PATH /opt/miniconda/bin:$PATH
ENV PATH /opt/velvet_1.2.10:$PATH

# Get basic ubuntu packages needed
RUN apt-get update -qq
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip

# Set up Miniconda environment for python2
RUN cd /opt;\
    wget http://repo.continuum.io/miniconda/Miniconda-3.3.0-Linux-x86_64.sh -O miniconda.sh;\
    chmod +x miniconda.sh;\
    ./miniconda.sh -p /opt/miniconda -b;\
    conda update --yes conda;\
    conda install --yes python=2.7

# Velvet for assembly
RUN apt-get install -qq zlib1g-dev
RUN cd /opt;\
    wget www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz -O velvet.tgz;\
    tar xf velvet.tgz;\
    cd velvet_1.2.10;\
    sed -i "s/MAXKMERLENGTH=31/MAXKMERLENGTH=128/" Makefile ;\
    make

# Bedtools2.17
RUN apt-get install -qq bedtools

# Picard tools 1.118
# To get fuse to work, I need the following (Issue here: https://github.com/dotcloud/docker/issues/514,
# solution here: https://gist.github.com/henrik-muehe/6155333).
ENV MRKDUP /opt/picard-tools-1.118/MarkDuplicates.jar
RUN apt-get install -qq libfuse2 openjdk-7-jre-headless
RUN cd /tmp ; apt-get download fuse
RUN cd /tmp ; dpkg-deb -x fuse_* .
RUN cd /tmp ; dpkg-deb -e fuse_*
RUN cd /tmp ; rm fuse_*.deb
RUN cd /tmp ; echo -en '#!/bin/bash\nexit 0\n' > DEBIAN/postinst
RUN cd /tmp ; dpkg-deb -b . /fuse.deb
RUN cd /tmp ; dpkg -i /fuse.deb
RUN cd /opt;\
    wget "http://downloads.sourceforge.net/project/picard/picard-tools/1.118/picard-tools-1.118.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fpicard%2Ffiles%2Fpicard-tools%2F1.118%2F&ts=1396879817&use_mirror=freefr" -O picard-tools-1.118.zip;\
    unzip picard-tools-1.118.zip

# Samtools 0.1.19
RUN apt-get install -qq samtools

# Bowtie2.1.0
RUN apt-get install -qq bowtie2

# Parallel 20130622-1
RUN apt-get install -qq parallel



# Install prodigal 2.60
RUN cd /opt;\
    wget --no-check-certificate https://prodigal.googlecode.com/files/Prodigal-2.60.tar.gz;\
    tar xf Prodigal-2.60.tar.gz;\
    cd Prodigal-2.60;\
    make;\
    ln -s /opt/Prodigal-2.60/prodigal /bin/prodigal

# Install R
RUN apt-get install -qq r-base

# Install R packages
RUN cd /opt;\
    RREPO='"http://cran.rstudio.com/"';\
    printf "install.packages(\"ggplot2\", repo=$RREPO)\ninstall.packages(\"reshape\",repo=$RREPO)\ninstall.packages(\"gplots\",repo=$RREPO)\ninstall.packages(\"ellipse\",repo=$RREPO)\ninstall.packages(\"grid\",repo=$RREPO)\ninstall.packages(\"getopt\",repo=$RREPO)" > dep.R;\
    Rscript dep.R

# Install python dependencies and fetch and install CONCOCT {{version}}
RUN cd /opt;\
    conda update --yes conda;\
    conda install --yes python=2.7 atlas cython numpy scipy biopython pandas pip scikit-learn pysam;\
    pip install bcbio-gff;\
    wget --no-check-certificate https://github.com/BinPro/CONCOCT/archive/{{version}}.tar.gz;\
    tar xf {{version}}.tar.gz;\
    cd CONCOCT-{{version}};\
    python setup.py install

ENV CONCOCT /opt/CONCOCT-{{version}}
ENV CONCOCT_TEST /opt/Data/CONCOCT-test-data
ENV CONCOCT_EXAMPLE /opt/Data/CONCOCT-complete-example

