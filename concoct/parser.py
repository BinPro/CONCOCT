import os
from argparse import ArgumentParser, ArgumentTypeError
from collections import namedtuple

import numpy as np

def set_random_state(seed):
    ERROR="'{0}' should be converatable to integer".format(seed)
    try:
        seed = int(seed)
        if seed < 0:
            raise ArgumentTypeError("'" + seed + "' should be >= 0")
        return seed
    except ValueError as e:
        raise ArgumentTypeError(ERROR)

def parse_split_pca(s):
    ERROR="'" + s + ("' is not a valid split pca proportion tuple. "
                     "Expected two positive integers <100")
    try:
        prop = s.split(",")
    except ValueError as e:
        raise ArgumentTypeError(ERROR)
    if not len(prop) == 2:
        raise ArgumentTypeError(ERROR)
    cov_prop = int(prop[0])
    comp_prop = int(prop([1]))
    return (cov_prop,comp_prop)

def arguments():
    parser = ArgumentParser()

    #Input files
    parser.add_argument('--coverage_file',
        help=("specify the coverage file, containing a table where each row "
              "correspond to a contig, and each column correspond to a sample. "
              "The values are the average coverage for this contig in that sample. "
              "All values are separated with tabs."))
    parser.add_argument('--composition_file',
        help=("specify the composition file, containing sequences in fasta format. "
              "It is named the composition file since it is used to calculate the "
              "kmer composition (the genomic signature) of each contig."))

    #Handle cluster number parsing
    parser.add_argument('-c', '--clusters', default=400, type=int,
      help='specify maximal number of clusters for VGMM, default 400.')
    #Kmer length, kmer count threshold and read length
    parser.add_argument('-k','--kmer_length', type=int, default=4,
        help='specify kmer length, defaults to tetramer')
    parser.add_argument('-l','--length_threshold', type=int, default=1000,
        help=("specify the sequence length threshold, contigs shorter than this "
              "value will not be used for fitting the model, but will be included "
              "in the final clustering results. Defaults to 1000"))
    parser.add_argument('-r','--read_length', type=int, default=100,
        help='specify read length for coverage, default 100')
    #Joined PCA or seperate PCA
    parser.add_argument('-s','--split_pca', action='store_true',
              help=('specify this flag to perform PCA for the composition '
               ' and coverage data seperately'))
    parser.add_argument('--coverage_percentage_pca', default=90, type=int,
                        help=('The percentatage of variance explained'
                              ' by the principal components for the'
                              ' coverage data, only considered if '
                              ' split pca is used'))
    parser.add_argument('--composition_percentage_pca', default=90, type=int,
                        help=('The percentatage of variance explained'
                              ' by the principal components for the'
                              ' composiiton data, only considered if '
                              ' split pca is used'))
    parser.add_argument('--total_percentage_pca', default=90, type=int,
                        help=('The percentage of variance explained'
                              ' by the principal components for the'
                              ' combined data, only considered if '
                              ' split pca is NOT used'))
    #Output
    parser.add_argument('-b', '--basename', default=os.curdir,
        help=("Specify the basename for files or directory where output"
              "will be placed. Path to existing directory or basename"
              "with a trailing '/' will be interpreted as a directory."
              "If not provided, current directory will be used."))
    parser.add_argument('-f','--force_seed',type=set_random_state, default=set_random_state(11),
                       help=('Specify an integer to use as seed for clustering. '
                             'You can specify 0 for random seed. The default seed '
                             'is 11.'))
    parser.add_argument('--no_cov_normalization', default=False, action="store_true",
                        help=("By default the coverage is normalized with regards to samples, "
                              "then normalized with regards of contigs and finally log transformed. "
                              "By setting this flag you skip the normalization and only do log "
                              "transorm of the coverage."))
    parser.add_argument('--no_total_coverage', default=False, action="store_true",
                        help=("By default, the total coverage is added as a new column in the coverage "
                              "data matrix, independently of coverage normalization but previous to "
                              "log transformation. Use this tag to escape this behaviour."))
    args  = parser.parse_args()
    # This can be changed to an or case when input of either case is supported individually
    if not (args.coverage_file and args.composition_file): 
        parser.error("No input data supplied, add file(s) using --coverage_file <cov_file> and/or "
                     "--composition_file <comp_file>")

    return args
