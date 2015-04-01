#!/usr/bin/env python
DESC= """
Transform a clustering result for contigs cutup with cut_up_fasta.py so that the original
contigs are assigned the majority vote among its parts.

prints result to stdout.

@author: alneberg
"""
import sys
import os
import argparse
import pandas as pd
from collections import Counter

def original_contig_name(s):
    """Transform s to the original contig name"""
    n = s.split(".")[-1]
    try:
        int(n)
    except:
        return s
    # Only small integers are likely to be 
    # indicating a cutup part.
    if int(n) < 1000:
        return ".".join(s.split(".")[:-1])
    else:
        # A large n indicates that the integer
        # was part of the original contig
        return s

def main(args):
    df = pd.read_table(args.clustering, sep=',', header=False, names=["contig_id", "cluster_id"])
    df['orig_contig_id'] = df.contig_id.apply(original_contig_name)
    majority_vote = {}
    for orig_contig_id, group_df in df.groupby('orig_contig_id'):
        c = Counter(group_df.cluster_id)
        majority_vote[orig_contig_id] = c.most_common(1)[0][0]
    output_s = pd.Series(majority_vote)
    output_s.to_csv(sys.stdout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument("clustering", help=("clustering_gtX.csv file from the CONCOCT output."
        "The contig names on the form contig_name.n, where n is an integer less than 1000 is "
        "assumed to originate from contig_name."))
    args = parser.parse_args()

    main(args)
