#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_true
from os.path import isdir
import os
import sys
import subprocess

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
        assert_true(isdir(tmp_dir_path),
                    msg = "Temporary directory not created")

