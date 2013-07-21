#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_true
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import re
import numpy as np

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'

CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """If tmp dir already exists, delete it"""
        if isdir(tmp_dir_path):
            self.tearDown()
        else:
            os.mkdir(tmp_dir_path)
        os.chdir(test_dir_path)

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            for f in os.listdir(d_path):
                f_path = os.path.join(d_path,f)
                os.remove(f_path)
            os.rmdir(d_path)

    def run_command(self,cov_file='coverage',comp_file='composition.fa'):
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                "CONCOCT test_data/{0} test_data/{1} '1','16' -o nose_tmp_output -c 3,5,1".format(cov_file,comp_file),
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def file_len(self,fh):
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_no_errors(self):
        self.run_command()
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        self.run_command()
        reg_o = re.compile('concoct_*')
        L = listdir(tmp_dir_path)
        tmp_dirs_before = filter(reg_o.match,L)
        assert_true(tmp_dirs_before!=0,
                    msg = "Temporary directory not created")
        assert_true(len(tmp_dirs_before)==1,
                    msg = "Other files starting with nose_tmp_output exists in test directory, please remove them and restart tests")

        # Rerun the concoct and see that new directory with unique
        # name is created
        self.run_command()

        assert_equal(self.c,0,
                     msg = "Error while running command a second time")
        L = listdir(tmp_dir_path)
        tmp_dirs_after = filter(reg_o.match,L)
        assert_true(len(tmp_dirs_after) > 1,
                    msg = "No unique output directory was created")

        assert_true(len(tmp_dirs_after) == 2,
                    msg = "Multiple output directories or files was created")

    def test_output_files_creation(self):
        self.run_command()
        for d in os.listdir(tmp_dir_path):
            d_p = os.path.join(tmp_dir_path,d)
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
                
    def test_threshold_functionality(self):
        self.run_command()
        for d in os.listdir(tmp_dir_path):
            d_p = os.path.join(tmp_dir_path,d)
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
            pca_m1 = np.loadtxt(pca_1)

        self.run_command(comp_file='composition_some_shortened.fa')
        for d in os.listdir(tmp_dir_path):
            d_temp = os.path.join(tmp_dir_path,d)
            if d_temp == d_p:
                continue
            d_p2 = d_temp
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
            pca_m2 = np.loadtxt(pca_2)

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
