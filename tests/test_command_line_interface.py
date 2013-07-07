#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal
import os
import sys
import fileinput

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestMultinomial(object):
    def setUp(self):
        pass
    def tearDown(self):
        """remove temporary output files"""
        if os.path.isdir('nose_tmp_output'):
            for f in os.listdir('nose_tmp_output'):
                os.remove(f)
            os.rmdir('nose_tmp_output')

    # testing function: log_probability
    def test_uniform_one_contig_prob(self):
        
        f = fileinput.input(os.path.join(data_path,"bambus2.scaffold.linear.fasta.one_contig"))
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        k = 4**4
        uniform_prob = np.ones((dna.DNA.kmer_hash_count))
        for i,cnt in dna_c.signature.items():
            uniform_prob[i] = 1./k
        log_prob = ml.log_probability(dna_c,uniform_prob)
        print log_prob
        assert_almost_equal(log_prob, -3791.05738056)
