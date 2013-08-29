# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 2013

@author: Johannes Alneberg
"""
import os
import sys
import pandas as p
import numpy as np
import logging

class Output(object):
    """
    Class to print out result information to their files
    """
    CONCOCT_PATH = None
    DT = None
    BIC_FILE = None
    ARGS_FILE = None
    PCA_FILE_BASE = None
    CLUSTERING_FILE_BASE = None
    MEANS_FILE_BASE = None
    VARIANCE_FILE_BASE = None
    RESPONSIBILITY_FILE_BASE = None
    @classmethod
    def __init__(self,basename,args):
        """
        Set output params and create output folders and bic.csv and args.txt
        """
        if os.path.isdir(basename):
            if basename[-1] == '/':
                self.CONCOCT_PATH = basename
            else:
                self.CONCOCT_PATH = basename + '/'
        elif basename[-1] == '/':
            basename_path = os.path.abspath(basename)
            os.mkdir(basename_path)
            self.CONCOCT_PATH = basename_path +'/'
        else:
            basename_path = os.path.abspath(basename)
            self.CONCOCT_PATH = basename+'_'

        print >> sys.stderr, "Results created at {0}".format(
            self.CONCOCT_PATH)
        self.BIC_FILE = self.CONCOCT_PATH + "bic.csv"
        self.ARGS_FILE = self.CONCOCT_PATH + "args.txt"
        self.ORIGINAL_FILE_BASE = self.CONCOCT_PATH + "original_data_gt{0}.csv"
        self.PCA_FILE_BASE = self.CONCOCT_PATH + \
            "PCA_transformed_data_gt{0}.csv"
        self.CLUSTERING_FILE_BASE = self.CONCOCT_PATH + "clustering{0}.csv"
        self.PCA_MEANS_FILE_BASE = self.CONCOCT_PATH + "pca_means{0}.csv"
        self.MEANS_FILE_BASE = self.CONCOCT_PATH + "means{0}.csv"
        self.VARIANCE_FILE_BASE = self.CONCOCT_PATH + "variance_gt{0}_dim{1}.csv"
        self.RESPONSIBILITY_FILE_BASE = self.CONCOCT_PATH + "responsibilities.csv"
        self.LOG_FILE_BASE = self.CONCOCT_PATH + 'concoct_log.txt'
        self.pipe = args.pipe

        logging.basicConfig(
            filename=self.LOG_FILE_BASE,
            level=logging.INFO,
            filemode='w', # Overwrites old log file
            format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
            )
        #Write header to bic.csv
        with open(self.BIC_FILE,"a+") as fh:
            print >> fh, "cluster_count,bic_value"
        with open(self.ARGS_FILE,"w+") as fh:
            print >> fh, args
    
    @classmethod
    def write_pca(self,transform,threshold,index):
        transform_df = p.DataFrame(transform,index=index)
        transform_df.to_csv(self.PCA_FILE_BASE.format(threshold))
    
    @classmethod
    def write_original_data(self,original,threshold):
        original.to_csv(self.ORIGINAL_FILE_BASE.format(threshold))

    
    @classmethod
    def write_clustering(self,dataframe,threshold_filter,threshold,c,pipe):
        if pipe:
            dataframe.clustering.to_csv(sys.stdout)
        dataframe.clustering.to_csv(
            self.CLUSTERING_FILE_BASE.format(""))
        dataframe[threshold_filter].clustering.to_csv(
            self.CLUSTERING_FILE_BASE.format(
                "_gt{0}".format(threshold)))
        
    @classmethod
    def write_bic(self,bics):
        bics.sort(key=lambda x: x[1])
        with open(self.BIC_FILE,"w+") as fh:
            for bic,c in bics:
                print >> fh, "{0},{1}".format(c,bic)
    
    @classmethod
    def write_cluster_pca_means(self,means,threshold,c):
        np.savetxt(
            self.PCA_MEANS_FILE_BASE.format("_gt{0}".format(threshold)),
            means,delimiter=',')

    @classmethod
    def write_cluster_means(self,means,threshold,c):
        np.savetxt(
            self.MEANS_FILE_BASE.format("_gt{0}".format(threshold)),
            means,delimiter=',')
            
    @classmethod
    def write_cluster_variance(self,var,threshold,i):
        np.savetxt(self.VARIANCE_FILE_BASE.format(threshold,i),
        var,delimiter=',')

    @classmethod
    def write_cluster_responsibilities(self,res,threshold,c):
        np.savetxt(self.RESPONSIBILITY_FILE_BASE,
        res,delimiter=',')
