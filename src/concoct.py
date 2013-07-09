#!/usr/bin/env python
from __future__ import division

import sys
import os
import re

import pandas as p
import pylab as pl
import matplotlib as mpl
import numpy as np

from itertools import tee, izip
from argparse import ArgumentParser,ArgumentTypeError
from sklearn.preprocessing import scale
from sklearn.mixture import GMM
from sklearn.decomposition import PCA
from datetime import datetime
from Bio import SeqIO


class Output(object):
    """
    Class to print out result information to their files
    """
    CONCOCT_PATH = None
    DT = None
    BIC_FILE = None
    ARGS_FILE = None
    
    @classmethod
    def __init__(self,args):
        """
        Set output parameters and create output folders and bic.csv and args.txt
        """
        self.DT = datetime.now().strftime("%y%m%d_%H%M")
        self.CONCOCT_PATH = os.path.join(args.o,"concoct_{0}".format(self.DT))
        os.makedirs(self.CONCOCT_PATH)
        print >> sys.stderr, "Results created in folder {0}".format(self.CONCOCT_PATH)
        self.BIC_FILE = os.path.join(self.CONCOCT_PATH,"bic.csv")
        self.ARGS_FILE = os.path.join(self.CONCOCT_PATH,"args.txt")
        #Write header to bic.csv
        with open(self.BIC_FILE,"a+") as fh:
            print >> fh, "cluster_count,bic_value"
        with open(self.ARGS_FILE,"w+") as fh:
            print >> fh, args
    
    @classmethod
    def write_cluster_method(self,c):
        print self.CONCOCT_PATH

        
def cluster(comp_file, cov_file, kmer_len, read_length, clusters_range, cov_range, split_pca, inits, iters, outdir, args):
    Output(args)
    #Composition
    #Generate kmer dictionary
    feature_mapping, nr_features = generate_feature_mapping(kmer_len)
    #Count lines in composition file
    count_re = re.compile("^>")
    seq_count = 0
    with open(comp_file) as fh:
        for line in fh:
            if re.match(count_re,line):
                seq_count += 1

    #Initialize with ones since we do pseudo count
    composition = np.ones(seq_count,nr_features)
    
    
    contigs_id = []
    for i,seq in enumerate(SeqIO.parse(comp_file,"fasta")):
        contigs_id.append(seq.id)
        for kmer_tuple in window(seq.seq.tostring().upper(),kmer_len):
            composition[i,feature_mapping["".join(kmer_tuple)]] += 1
    composition = p.DataFrame(composition,index=contigs_id,dtype=float)
    threshold_filter = composition.sum(axis=1) > 1000
    #log(p_ij) = log((X_ij +1) / rowSum(X_ij+1))
    composition = np.log(composition.divide(composition.sum(axis=1),axis=0))
    
    #Coverage import, file has header and contig ids as index
    cov = p.read_table(cov_file,header=0,index_col=0)
        
    
    #cov = scale(cov.ix[:,cov_range[0]:cov_range[1]])
    cv_type='full'
    
    
    for c in clusters_range:
        # Fit a mixture of gaussians with EM
        gmm = GMM(n_components=c, covariance_type=cv_type)
        gmm.fit(X)
        labels = gmm.predict(X)
        bic = gmm.bic(X)
        convergence = gmm.converged_
        sys.stderr.write('Convergence: ' + str(convergence) +'\n')
        class_series = p.Series(labels,index=df.index)
        setting = str(c)+"_" + str(cv_type)
        
        class_series.to_csv(outdir+"_" +setting+str(gmm.bic(X))+ '.csv')

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in xrange(i):
            next(el, None)
    return izip(*els)

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        kmer = ''.join(kmer)
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = ''.join([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter+1

def parse_cluster_list(cc_string):
    ERROR="'" + cc_string + "' is not a valid range of number. Expected forms like '20,100,2'."
    try:
        cc=map(int,cc_string.split(","))
    except ValueError as e:
        raise ArgumentTypeError(ERROR)    
    if not len(cc) == 3:
        raise ArgumentTypeError(ERROR)
    first,last,step = cc
    return xrange(first, last+1, step)
    
def parse_coverage_columns(cov_string):
    ERROR="'" + cov_string + "' is not valid. Expected 'first_column_name,last_column_name'."
    try:
        cov = cov_string.split(",")
    except ValueError as e:
        raise ArgumentTypeError(ERROR)
    if not len(cov) == 2:
        raise ArgumentTypeError(ERROR)
    return cov

def parse_taxonomy_cluster_list(tax_file):
    raise NotImplementedError("This functionality has not been added yet. Please use -c and specify range")

def arguments():
    parser = ArgumentParser()
    #Handle cluster number parsing
    cluster_count = parser.add_mutually_exclusive_group()
    cluster_count.add_argument('-c', nargs="+", default=range(20,101,2), type=parse_cluster_list,
        help='specify range of clusters to try out on format first,last,step. default 20,100,2.')
    cluster_count.add_argument('-t', type=parse_taxonomy_cluster_list,
        help='specify a taxonomy file to estimate species number from (X). Will use range X*0.5,X*1.5,2')
        
    parser.add_argument('coverage_file',
        help='specify the coverage file')
    parser.add_argument('composition_file',
        help='specify the composition file')
    parser.add_argument('-k', type=int, default=4,
        help='specify kmer length, defaults to tetramer')
    parser.add_argument('-r', type=int, default=100,
        help='specify read length for coverage, default 100')
    parser.add_argument('-n', type=parse_coverage_columns, default=None,
        help='specify the first and last column names for continuous coverage range of read counts as first,last')
    parser.add_argument('-s', type=bool, default=False, action="store_true",
        help='specify this flag to first do PCA for the composition and using that component number\
              that explaines 90% of variance for the coverage as well. Default join composition and\
              coverage before PCA.')
    
    parser.add_argument('-e', type=int, default=5,
        help='How often to initialize each cluster count. default 5 times')
    parser.add_argument('-i', type=int, default=100,
        help='Maximum number of iterations if convergance not achieved')

    parser.add_argument('-o', default=os.curdir,
        help='specify the output directory, if not provided current directory used. All files will\
              be created in the folder/CONCOCT_YYMMDD_HHMM folder')
    
    
    return parser.parse_args()

        
if __name__=="__main__":
    args = arguments()
    results = cluster(args.composition_file, args.coverage_file, args.k, args.r, args.c, args.n, \
                      args.s, args.e, args.i, args.o, args)
    
