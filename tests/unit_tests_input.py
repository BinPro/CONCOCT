#!/usr/bin/env python
from nose.tools import assert_equal, assert_true
import numpy as np
import pandas as p
from concoct.input import _normalize_per_sample, _normalize_per_contig

class TestInput(object):
    def setUp(self):
        self.C = p.DataFrame(np.array([[0., 0.7], [5.5, .7]]))

    def test_normalize_per_contig(self):
        C_norm = _normalize_per_contig(self.C)
        
        C_correct = p.DataFrame(np.array([[0., 1.],[5.5/6.2, 0.7/6.2]]))
        assert_true(np.linalg.norm(C_norm-C_correct) < 0.0001)

    def test_normalize_per_samples(self):
        C_norm = _normalize_per_sample(self.C)

        C_correct = p.DataFrame(np.array([[0., 0.5],[1,0.5]]))
        assert_true(np.linalg.norm(C_norm-C_correct) < 0.0001)

