import sys
import os
from os.path import join as ospj
from nose.tools import assert_equal
import pandas as pd

import concoct.utils.dir_utils as dir_utils

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data", "scg_bins"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, 'nose_tmp_output')
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, 'extract_scg_bins')
SCRIPT_PATH = ospj(TEST_DIR_PATH, '..')

# Add script dir to python path to import functions
sys.path.append(SCRIPT_PATH)
from extract_scg_bins import get_approved_bins, sum_bases_in_bins, \
    get_winning_bins

CWD = os.getcwd()


class TestDnaDiff(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        dir_utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temporary output files"""
        dir_utils.rm_rf(TMP_DIR_PATH)

    def test_get_approved_bins(self):
        """Test get_approved_bins"""
        df = get_approved_bins(ospj(DATA_PATH, "sample1_gt500_scg.tsv"),
                max_missing_scg=2, max_multicopy_scg=4)
        assert_equal(6, int(df.Cluster))

    def test_sum_bases_in_bins(self):
        """Test sum_bases_in_bins"""
        scg_tsv = ospj(DATA_PATH, "sample1_gt500_scg.tsv")
        b = sum_bases_in_bins(pd.read_csv(scg_tsv, sep="\t"),
                ospj(DATA_PATH, "sample1_gt500.fa"))
        assert_equal(16170706, b)
        df = get_approved_bins(ospj(DATA_PATH, "sample1_gt500_scg.tsv"),
                max_missing_scg=2, max_multicopy_scg=4)
        b = sum_bases_in_bins(df, ospj(DATA_PATH, "sample1_gt500.fa"))
        assert_equal(5432170, b)

    def test_get_winning_bins(self):
        """Test get_winning_bins"""
        scg_tsvs = [ospj(DATA_PATH, p) for p in ["sample1_gt300_scg.tsv",
            "sample1_gt500_scg.tsv"]]
        fasta_files = [ospj(DATA_PATH, p) for p in ["sample1_gt300.fa",
            "sample1_gt500.fa"]]
        winning_index, df = get_winning_bins(scg_tsvs, fasta_files,
                max_missing_scg=2, max_multicopy_scg=4)
        assert_equal(1, winning_index)
        winning_index, df = get_winning_bins(list(reversed(scg_tsvs)),
            list(reversed(fasta_files)), max_missing_scg=2,
            max_multicopy_scg=4)
        assert_equal(0, winning_index)
