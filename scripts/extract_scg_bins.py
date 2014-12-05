#!/usr/bin/env python
"""
Extract bins with given SCG criteria.
"""
import pandas as pd

def get_approved_bins(cog_tsv, max_missing_scg=2, max_multicopy_scg=4):
    cogdf = pd.read_csv(cog_tsv, sep="\t")
    # number of multicopy genes
    multi_scgs = cogdf.iloc[:,range(3, len(cogdf.columns))].apply(lambda x: x > 1).sum(axis=1)
    # number of missing scgs
    miss_scgs = cogdf.iloc[:,range(3, len(cogdf.columns))].apply(lambda x: x == 0).sum(axis=1)
    app_bins = cogdf[(miss_scgs <= max_missing_scg) & (multi_scgs <= max_multicopy_scg)].loc[:,["Cluster", "Contigs"]]

    return app_bins


def main(clustering_csv, cog_tsv, fasta_file):


if __name__ == "__main__":
    main(*parse_input())
