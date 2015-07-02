#!/usr/bin/env python
"""A script to create a concoct input table from kallisto abundance.txt output files"""

import argparse
import pandas as pd
import os
import sys

def main():
    sample_dfs = []
    for i, sample in enumerate(args.quantfiles):
        if args.samplenames:
            samplename = args.samplenames[i]
        else:
            samplename = os.path.basename(sample)        
        sample_df = pd.read_table(sample, index_col=0)
        
        sample_dfs.append((samplename, sample_df))
    kallisto_df = pd.DataFrame(index=sample_df.index)

    for sample, sample_df in sample_dfs:
        kallisto_df['kallisto_coverage_{0}'.format(sample)] = 200*sample_df['est_counts'].divide(sample_df['length'])
    kallisto_df.to_csv(sys.stdout, sep="\t", float_format="%.6f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quantfiles", nargs='+', help="Kallisto abundance.txt files")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each, Should be the same order and the same number as the abundance.txt files")
    args = parser.parse_args()

