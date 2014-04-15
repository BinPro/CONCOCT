# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v0.2.2
# VERSION 0.2.2
# 
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v0.2.2 installed.
#
# To login to the docker with a shared directory from the host do:
#
# sudo docker run -v /my/host/shared/directory:/my/docker/location -i -t binnisb/concoct_0.2.2 /bin/bash
#
# This environment does not set up the assembler and preprocessing
# for concoct. We will be creating another docker for that.

FROM ubuntu:13.10
MAINTAINER CONCOCT developer group, concoct-support@lists.sourceforge.net

ENV PATH /opt/miniconda/bin:$PATH

# Get basic ubuntu packages needed
RUN apt-get update -qq;\
    apt-get install -qq wget build-essential libgsl0-dev 

# Set up Miniconda environment for python2
RUN cd /opt;\
    wget http://repo.continuum.io/miniconda/Miniconda-3.3.0-Linux-x86_64.sh -O miniconda.sh;\
    chmod +x miniconda.sh;\
    ./miniconda.sh -p /opt/miniconda -b

# Install python dependencies and fetch and install CONCOCT 0.2.2
RUN cd /opt;\
    conda update --yes conda;\
    conda install --yes python=2.7  atlas cython numpy scipy biopython pandas pip scikit-learn;\
    wget --no-check-certificate https://github.com/BinPro/CONCOCT/archive/0.2.2.tar.gz;\
    tar xf 0.2.2.tar.gz;\
    cd CONCOCT-0.2.2;\
    python setup.py install



