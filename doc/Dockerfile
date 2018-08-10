# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v0.5.0
# VERSION 0.5.0
#
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v0.5.0 installed.
#
# To login to the docker with a shared directory from the host do:
#
# sudo docker run -v /my/host/shared/directory:/my/docker/location -i -t binnisb/concoct_0.5.0 /bin/bash
#

FROM continuumio/miniconda
MAINTAINER CONCOCT developer group, concoct-support@lists.sourceforge.net
COPY . /opt/CONCOCT

# Get basic ubuntu packages needed
RUN apt-get update -qq
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip bedtools

RUN conda update --yes conda;\
    conda install --yes python=2.7;\
    pip install --upgrade pip

# Install python dependencies and fetch and install CONCOCT 0.5.0
RUN cd /opt;\
    conda install --yes python=2.7 atlas cython numpy scipy biopython pandas pip scikit-learn pysam tqdm;\
    pip install bcbio-gff
    
#    wget --no-check-certificate https://github.com/BinPro/CONCOCT/archive/0.5.0.tar.gz;\
#    tar xf 0.5.0.tar.gz;\
#    cd CONCOCT-0.5.0;\
#    python setup.py install

RUN cd /opt/CONCOCT/;\
    python setup.py install

RUN cd /opt/CONCOCT/;\
    nosetests
