#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_true
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import re

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'

CWD = os.getcwd()

CONCOCT_CALL = """
CONCOCT test_data/coverage -o nose_tmp_output -c 3
"""
class TestCMD(object):
    def setUp(self):
        """If tmp dir already exists, delete it"""
        self.tearDown()
        os.chdir(test_dir_path)
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                CONCOCT_CALL,shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode
    def tearDown(self):
        """remove temporary output files"""
        if isdir(tmp_dir_path):
            for f in os.listdir(tmp_dir_path):
                os.remove(f)
            os.rmdir(tmp_dir_path)

    def test_no_errors(self):
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        reg_o = re.compile(tmp_dir_path)
        L = listdir(test_dir_path)
        tmp_dirs_before = filter(reg_o.match,L)

        assert_true(isdir(tmp_dir_path),
                    msg = "Temporary directory not created")
        assert_true(len(tmp_dirs_before)==1,
                    msg = "Other files starting with nose_tmp_output exists in test directory, please remove them and restart tests")

        # Rerun the concoct and see that new directory with unique
        # name is created
        self.op = subprocess.check_output(
            CONCOCT_CALL,shell=True)
        L = listdir(test_dir_path)
        tmp_dirs_after = filter(reg_o.match,L)
        assert_true(len(tmp_dirs_after) > 1,
                    msg = "No unique output directory was created")

        assert_true(len(tmp_dirs_after) == 2,
                    msg = "Multiple output directories or files was created")

    def test_output_files_creation(self):
        assert_true(
            isfile(tmp_dir_path + '/bic'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/large_contigs_clustering.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/large_contigs_cluster_means.csv'),
            msg='Large contigs cluster means file is not created'
            )

        assert_true(
            isfile(tmp_dir_path+ '/large_contigs_cluster_variance.csv'),
            msg='Large contigs cluster variance file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/large_contigs_clustering.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/large_contigs_responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/pca_data.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(tmp_dir_path+ '/original_data.csv'),
            msg='Original data file is not created'
            )
                
