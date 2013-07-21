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
CONCOCT test_data/coverage test_data/composition.fa '1','16' -o nose_tmp_output -c 3,5,1
"""
class TestCMD(object):
    def setUp(self):
        """If tmp dir already exists, delete it"""
        if isdir(tmp_dir_path):
            self.tearDown()
        else:
            os.mkdir(tmp_dir_path)
        os.chdir(test_dir_path)
        self.run_command()

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            for f in os.listdir(d_path):
                f_path = os.path.join(d_path,f)
                os.remove(f_path)
            os.rmdir(d_path)

    def run_command(self):
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                CONCOCT_CALL,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def test_no_errors(self):
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        reg_o = re.compile('concoct_*')
        L = listdir(tmp_dir_path)
        tmp_dirs_before = filter(reg_o.match,L)
        sys.stderr.write(str(tmp_dirs_before) + '\n')
        assert_true(tmp_dirs_before!=0,
                    msg = "Temporary directory not created")
        assert_true(len(tmp_dirs_before)==1,
                    msg = "Other files starting with nose_tmp_output exists in test directory, please remove them and restart tests")

        # Rerun the concoct and see that new directory with unique
        # name is created
        self.run_command()

        assert_equal(self.c,0,
                     msg = "Error while running command a second time")
        L = listdir(test_dir_path)
        tmp_dirs_after = filter(reg_o.match,L)
        assert_true(len(tmp_dirs_after) > 1,
                    msg = "No unique output directory was created")

        assert_true(len(tmp_dirs_after) == 2,
                    msg = "Multiple output directories or files was created")

    def test_output_files_creation(self):
        for d in os.listdir(tmp_dir_path):
            d_p = os.path.join(tmp_dir_path,d)
            sys.stderr.write(d_p+'\n')
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
                
