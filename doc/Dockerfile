# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v0.2.1
# VERSION 0.2.1
# 
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v0.2.1 installed.
#
# To login to the docker with a shared directory from the host do:
#
# sudo docker run -v /my/host/shared/directory:/my/docker/location -i -t binnisb/concoct_0.2.1 /bin/bash
#
# This environment does not set up the assembler and preprocessing
# for concoct. We will be creating another docker for that.


FROM ubuntu:13.10
MAINTAINER CONCOCT developer group, concoct-support@lists.sourceforge.net

RUN apt-get update -qq;\
    apt-get install -qq wget git build-essential libgsl0-dev gsl-bin

RUN cd /opt;\
    wget http://repo.continuum.io/miniconda/Miniconda-3.3.0-Linux-x86_64.sh -O miniconda.sh;\
    chmod +x miniconda.sh;\
    ./miniconda.sh -p /opt/miniconda -b

ENV PATH /opt/miniconda/bin:$PATH

RUN cd /opt;\
    conda update --yes conda;\
    conda install --yes python=2.7  atlas cython numpy scipy biopython pandas pip scikit-learn;\
    git clone https://github.com/BinPro/CONCOCT.git;\
    cd CONCOCT;\
    git checkout tags/0.2.1;\
    python setup.py install





