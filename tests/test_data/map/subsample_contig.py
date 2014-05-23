#!/usr/bin/env python
"""extract_fasta_bins.py 

Extract a fasta file for each cluster from a concoct result file.
"""

import argparse
import sys
import os
import numpy as np
from Bio import SeqIO
from Bio import Seq

def main(args):
    all_seqs = []
    for i, seq in enumerate(SeqIO.parse(args.fasta_file, "fasta")):
        for n in range(args.number_of_seqs):
            rand_i = np.random.randint(0, len(seq) - args.subseq_length + 1)
            subseq = seq[rand_i:rand_i + args.subseq_length] 
            subseq.id += "_read_{0}".format(n)
            all_seqs.append(subseq)
    SeqIO.write(all_seqs, sys.stdout, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_file", help="Input Fasta file.")
    parser.add_argument("number_of_seqs", type=int, help="Number of subsequences per input sequence")
    parser.add_argument("subseq_length", type=int, help="Number of bases in each subsequence")
    args = parser.parse_args()

    main(args)
