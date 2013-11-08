import os
from argparse import ArgumentParser, ArgumentTypeError
from collections import namedtuple

import numpy as np
import multiprocessing

def set_random_state(seed):
    ERROR="'{0}' should be converatable to integer".format(seed)
    try:
        seed = int(seed)
        if seed < 0:
            raise ArgumentTypeError("'" + seed + "' should be >= 0")
        return seed
    except ValueError as e:
        raise ArgumentTypeError(ERROR)

def get_max_n_processors(n_procs):
    #-------------------------------------------------------------------------------
    # MPI setup
    #-------------------------------------------------------------------------------
    MpiParams = namedtuple("MpiParams","comm use_mpi size rank")
    if n_procs:
        try:
            n_procs = int(n_procs)
        except ValueError:
            raise ArgumentTypeError("{0} should be convertable to integer".format(n_procs))
    try:
        if not os.environ.has_key('OMPI_COMM_WORLD_RANK'):
            raise ImportError
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        use_mpi = True
        size =  n_procs if n_procs else comm.Get_size() 
        rank = comm.Get_rank()
    except ImportError:
        comm = None
        use_mpi = False
        rank = None
        size = n_procs if n_procs else multiprocessing.cpu_count()

    return MpiParams(comm,use_mpi,size,rank)

def parse_cluster_list(cc_string):
    ERROR="'" + cc_string + ("' is not a valid range of number. Expected "
                             "forms like '20,100,2'.")
    try:
        first, last, step = map(int,cc_string.split(","))
    except ValueError as e:
        raise ArgumentTypeError(ERROR)
    except Exception as e:
        raise ArgumentTypeError(ERROR)
    return xrange(first, last+1, step)


def parse_taxonomy_cluster_list(tax_file):
    raise NotImplementedError(("This functionality has not been added yet. "
                               "Please use -c and specify range"))


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
        help='specify the coverage file')
    parser.add_argument('--composition_file',
        help='specify the composition file')

    #Handle cluster number parsing
    cluster_count = parser.add_mutually_exclusive_group()
    cluster_count.add_argument('-c', '--clusters', default=range(20,101,2), 
                               type=parse_cluster_list,
                               help=('specify range of clusters to try out'
                                     ' on format first,last,step.'
                                     ' default 20,100,2.'))
    #cluster_count.add_argument('-t', type=parse_taxonomy_cluster_list,
    #help='specify a taxonomy file to estimate species number from (X). \
    #      Will use range X*0.5,X*1.5,2')



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
    #Clustering Parameters
    parser.add_argument('-e', '--executions',type=int, default=5,
        help='How often to initialize each cluster count. default 5 times')
    parser.add_argument('-i', '--iterations',type=int, default=1000,
        help='Maximum number of iterations if convergance not achieved')
    #Output
    parser.add_argument('-b', '--basename', default=os.curdir,
        help=("Specify the basename for files or directory where output"
              "will be placed. Path to existing directory or basename"
              "with a trailing '/' will be interpreted as a directory."
              "If not provided, current directory will be used."))
    parser.add_argument('-p', '--pipe', default=False, action="store_true",
                        help=('Add this tag if the main result file should be'
                              'printed to stdout. Useful for pipeline use'))
    parser.add_argument('-m','--max_n_processors',type=get_max_n_processors, 
                        default=get_max_n_processors(None),
                        help=('Specify the maximum number of processors CONCOCT is allowed to use, '
                            'if absent, all present processors will be used. Default is to  '
                            'assume MPI execution, allowing execution over multiple nodes, '
                            'with a fallback to standard python multiprocessing on a single '
                            'machine using available cores. It is recommended to install mpi '
                            'and mpi4py if execution over multiple '
                            'nodes is required. To run with MPI use call it with mpirun -np N '
                            'or equivalent'))
    parser.add_argument('-f','--force_seed',type=set_random_state, default=set_random_state(11),
                       help=('Specify an integer to use as seed for clustering. '
                             'You can specify 0 for random seed. The default seed '
                             'is 11.'))
    parser.add_argument('--covariance_type', default="full", 
                        choices=['full','diag'], 
                        help=("Choose the shape of the covariance matrix for "
                              "the GMM:s used in clustering."))

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
