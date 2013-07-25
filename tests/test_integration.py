#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_true
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import re
import numpy as np
import pandas as p


file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'
tmp_basename_dir = tmp_dir_path + '/1'
tmp_basename_dir2 = tmp_dir_path + '/2'
tmp_basename_file = tmp_dir_path + '/file'


CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """If tmp dir already exists, delete it"""
        self.tearDown()
        os.mkdir(tmp_basename_dir)
        os.chdir(test_dir_path)

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            try:
                os.remove(d_path)
            except:
                for f in os.listdir(d_path):
                    f_path = os.path.join(d_path,f)
                    os.remove(f_path)
                os.rmdir(d_path)
        assert os.listdir(tmp_dir_path) == []


    def run_command(self,cov_file='coverage',comp_file='composition.fa',tags=[],basename='nose_tmp_output/1'):
        call_string = "CONCOCT test_data/{0} test_data/{1} --basename {2} -c 3,5,1".format(cov_file,comp_file,basename)
        for tag in tags:
            call_string += " " + tag
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def file_len(self,fh):
        i=0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def md5sum(self,fh):
        infile = open("filename", 'rb')
        content = infile.read()
        infile.close()
        m = hashlib.md5() 
        m.update(content)
        return m.hexdigest()
 
    def test_no_errors(self):
        self.run_command()
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        self.run_command()
        assert_true(isdir(tmp_basename_dir),
                    msg = "Temporary directory not created")
        m_time_first = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        
        # Rerun the concoct and see that the directory is overwritten
        self.run_command()
        m_time_second = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        assert_true(m_time_first != m_time_second,
                     msg = "basename dir is not overwritten")
        L = listdir(tmp_dir_path)
        assert_true(len(L) == 1,
                    msg = "Multiple output directories or files was created")

        # File creation
        self.run_command(basename=tmp_basename_file)
        assert_true(isfile(tmp_basename_file+'_clustering.csv'),
                    msg = "Clustering file is not created, when file is used as basename")
        L = listdir(tmp_dir_path)
        assert_true(len(L) == 12,
                    msg = "Wrong number of output files")

    def test_output_files_creation(self):
        # dir as basename
        self.run_command()
        d_p = os.path.join(tmp_basename_dir)
        assert_true(
            isfile(d_p+ '/bic.csv'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/means_gt1000.csv'),
            msg='Large contigs cluster means file is not created'
            )
        
        assert_true(
            isfile(d_p+ '/variance_gt1000_dim1.csv'),
            msg='Large contigs cluster variance file is not created'
            )
        assert_true(
            isfile(d_p+ '/responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(d_p+ '/original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        
        # dir as file
        self.run_command(basename=tmp_basename_file)
        d_p = tmp_basename_file +'_'
        assert_true(
            isfile(d_p +'bic.csv'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(d_p+ 'clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(d_p+ 'means_gt1000.csv'),
            msg='Large contigs cluster means file is not created'
            )

        assert_true(
            isfile(d_p+ 'variance_gt1000_dim1.csv'),
            msg='Large contigs cluster variance file is not created'
            )
        assert_true(
            isfile(d_p+ 'responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(d_p+ 'clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(d_p+ 'PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(d_p+ 'original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
                
    def test_threshold_functionality(self):
        self.run_command()
        d_p = tmp_basename_dir
        od_1 = d_p+'/original_data_gt1000.csv'
        pca_1 = d_p+'/PCA_transformed_data_gt1000.csv'
        var_1 = d_p+'/variance_gt1000_dim1.csv'
        means_1 = d_p+'/means_gt1000.csv'
        clust_gt_1 = d_p+'/clustering_gt1000.csv'
        clust_1 = d_p+'/clustering.csv'
        odl_1 = self.file_len(od_1)
        varl_1= self.file_len(var_1)
        meansl_1= self.file_len(means_1)
        clust_gtl_1= self.file_len(clust_gt_1)
        clustl_1 = self.file_len(clust_1)
        
        pca_df1 = p.io.parsers.read_table(pca_1)
        pca_m1 = pca_df1.to_records()

        self.run_command(comp_file='composition_some_shortened.fa',
                         basename=tmp_basename_dir2+'/')
        d_p2 = tmp_basename_dir2
        od_2 = d_p2+'/original_data_gt1000.csv'
        pca_2 = d_p2+'/PCA_transformed_data_gt1000.csv'
        var_2 = d_p2+'/variance_gt1000_dim1.csv'
        means_2 = d_p2+'/means_gt1000.csv'
        clust_gt_2 = d_p2+'/clustering_gt1000.csv'
        clust_2 = d_p2+'/clustering.csv'
        odl_2 = self.file_len(od_2)
        varl_2= self.file_len(var_2)
        meansl_2= self.file_len(means_2)
        clust_gtl_2= self.file_len(clust_gt_2)
        clustl_2 = self.file_len(clust_2)
        pca_df2 = p.io.parsers.read_table(pca_2)
        pca_m2 = pca_df2.to_records()
        
        assert_true(odl_1!=odl_2,
                    msg='Original data have the same lengths')
        assert_true(varl_1==varl_2,
                    msg='Variance files does not have the same lengths')
        assert_true(meansl_1==meansl_2,
                    msg='Means files does not have the same lengths')
        assert_true(clust_gtl_1!=clust_gtl_2,
                    msg='Filtered clustering files have the same lengths')
        assert_true(pca_m1.shape!=pca_m2.shape,
                    msg='PCA transformed data has the same shapes')
        assert_true(clustl_1==clustl_2,
                    msg='Clustering files does not have the same lengths')
        assert_true(clust_gtl_2!=clustl_2,
                    msg='Filtered clustering file and full have the same lengths')

    def test_piping_functionality(self):
        f1 = tmp_basename_dir+'/stdout_capture'
        self.run_command(tags=['--pipe','> {0}/stdout_capture'.format(tmp_basename_dir)])

        f1 = open(tmp_basename_dir+'/stdout_capture','rb')
        f1_content = f1.read()
        f1.close()
        f2 = open(tmp_basename_dir + '/clustering.csv','rb')
        f2_content = f2.read()
        f2.close()
        import filecmp
        assert_true(len(f1_content)==len(f2_content),
                    msg='stdout and clustering file is not equal')

    def test_bic_sorted(self):
        self.run_command()
        bic = p.io.parsers.read_table(tmp_basename_dir+'/bic.csv',sep=',',index_col=0,header=None)
        assert_true(max(bic.index) == 5,
                    msg='BIC columns are probably mixed up')
        index_l = list(bic.index)
        assert_true(index_l==[3,4,5],
                    msg='BIC file is not sorted')
        # Run command again, to see that file is overwritten
        self.run_command()
        bic = p.io.parsers.read_table(tmp_basename_dir+'/bic.csv',sep=',',index_col=0,header=None)
        assert_true(len(bic)==3,
                    'BIC file is probably appended to')
