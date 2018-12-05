#!/usr/bin/env python
"""
With contigs cutup with cut_up_fasta.py as input, sees to that the consequtive
parts of the original contigs are merged.

prints result to stdout.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
from collections import defaultdict, Counter

def original_contig_name_special(s):
    """Transform s to the original contig name according to the special Baltic Mag paper"""
    n = s.split(".")[-1].split('_')[0]
    try:
        int(n)
    except:
        return s, 0, None
    # Only small integers are likely to be 
    # indicating a cutup part.
    if int(n) < 1000:

        return ".".join(s.split(".")[:-1]), int(n), "_".join(s.split(".")[-1].split('_')[1:])
    else:
        # A large n indicates that the integer
        # was part of the original contig
        return s, 0, None

def main(args):
    all_seqs = {}
    all_originals = defaultdict(dict)
    with open(args.cutup_clustering_result, 'r') as ifh:
        for line in ifh:
            line = line.strip()
            contig_id, cluster_id = line.split(',')
            original_contig_name, part_id, left_over = original_contig_name_special(contig_id)
        
            all_originals[original_contig_name][part_id] = cluster_id

    merged_contigs_stack = []
    
    for original_contig_id, part_ids_d in all_originals.iteritems():
        if len(part_ids_d) > 1:
            c = Counter(part_ids_d.values())
            cluster_id = c.most_common(1)[0][0]
            c_string = [(a,b) for a, b in c.iteritems()]
            if len(c.values()) > 1:
                sys.stderr.write("{}\t{}, chosen: {}\n".format(original_contig_id, c_string, cluster_id))
            else:
                sys.stderr.write("{}\t{}\n".format(original_contig_id, c_string))
        else:
            cluster_id = part_ids_d.values()[0]     

        sys.stdout.write("{},{}\n".format(original_contig_id, cluster_id))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cutup_clustering_result", help=("Input cutup clustering result."))
    args = parser.parse_args()

    main(args)
