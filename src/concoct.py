#!/usr/bin/env python
import fileinput
import sys
import os
import re
import itertools

import pandas as p
import sklearn.mixture as mixture
import pylab as pl
import matplotlib as mpl
import numpy as np

from argparse import ArgumentParser,ArgumentTypeError
from scipy import linalg
from sklearn import preprocessing



def main(coverage_file,output_base,n_components,header=None):
    df = p.read_csv(coverage_file,header=header,index_col=0)
    X = df.values
    X = preprocessing.scale(X)
    cv_type='full'

    # Fit a mixture of gaussians with EM
    gmm = mixture.GMM(n_components=n_components, covariance_type=cv_type)
    gmm.fit(X)
    labels = gmm.predict(X)
    bic = gmm.bic(X)
    convergence = gmm.converged_
    sys.stderr.write('Convergence: ' + str(convergence) +'\n')
    class_series = p.Series(labels,index=df.index)
    setting = str(n_components)+"_" + str(cv_type)
    
    class_series.to_csv(output_base+"_" +setting+str(gmm.bic(X))+ '.csv')
    
    

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

def parse_taxonomy_cluster_list(tax_file):
    raise NotImplementedError("This functionality has not been added yet. Please use -c and specify range")

def arguments():
    parser = ArgumentParser()
    parser.add_argument('coverage_file',
        help='specify the coverage file')
    parser.add_argument('composition_file',
        help='specify the composition file')
    parser.add_argument('-k', type=int, default=4,
        help='specify kmer length, defaults to tetramer')

    parser.add_argument('-r', type=int, default=5,
        help='How often to initialize each cluster count. default 5 times')
    parser.add_argument('-i', type=int, default=100,
        help='Maximum number of iterations if convergance not achieved')

    parser.add_argument('-o', 
        help='specify the output directory, if not provided current directory used. All files will\
              be created in the folder/CONCOCT_YYMMDD_HHMM folder')
    
    #Handle cluster number parsing
    cluster_count = parser.add_mutually_exclusive_group()
    cluster_count.add_argument('-c', type=int,nargs="+", default=xrange(20,101,2), type=parse_cluster_list,
        help='specify range of clusters to try out on format first,last,step. default 20,100,2.')
    cluster_count.add_argument('-t', type=parse_taxonomy_cluster_list,
        help='specify a taxonomy file to estimate species number from (X). Will use range X*50%,X*150%,2')
    
    return parser.parse_args()
    
if __name__=="__main__":
    args = arguments()
    print args

