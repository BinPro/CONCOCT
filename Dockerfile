# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v1.0.0
# VERSION 1.0.0
#
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v1.0.0 installed.
#
# To login to the docker with a shared directory from the host do:
#
# docker run -v /my/host/shared/directory:/my/docker/location -i -t alneberg/concoct_1.0.0 /bin/bash
#

FROM ubuntu:18.04
COPY . /opt/CONCOCT

# Get basic ubuntu packages needed
RUN apt-get update -qq
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip bedtools python-pip

RUN pip install --upgrade pip

# Install python dependencies and fetch and install CONCOCT 1.0.0
RUN cd /opt/CONCOCT;\
    pip install -r requirements.txt;\
    
#    wget --no-check-certificate https://github.com/BinPro/CONCOCT/archive/1.0.0.tar.gz;\
#    tar xf 1.0.0.tar.gz;\
#    cd CONCOCT-1.0.0;\
#    python setup.py install

RUN cd /opt/CONCOCT/;\
    python setup.py install

RUN cd /opt/CONCOCT/;\
    nosetests
