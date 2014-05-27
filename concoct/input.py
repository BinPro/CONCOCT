import re
import logging

import numpy as np
import pandas as p

from itertools import product, tee, izip

from Bio import SeqIO

def load_data(args):
    composition, contig_lengths, threshold_filter = load_composition(
        args.composition_file, 
        args.kmer_length, 
        args.length_threshold
        )

    if args.coverage_file:
        cov, cov_range = load_coverage(
            args.coverage_file,
            contig_lengths,
            args.no_cov_normalization,
            add_total_coverage = (not args.no_total_coverage),
            read_length = args.read_length
            )
    else:
        cov, cov_range = None, None

    return composition, contig_lengths, threshold_filter, cov, cov_range

def load_composition(comp_file,kmer_len,threshold):
    #Composition
    #Generate kmer dictionary
    feature_mapping, nr_features = generate_feature_mapping(kmer_len)
    #Count lines in composition file
    count_re = re.compile("^>")
    seq_count = 0
    with open(comp_file) as fh:
        for line in fh:
            if re.match(count_re,line):
                seq_count += 1

    # Initialize with ones since we do pseudo count, we have i contigs as rows
    # and j features as columns
    composition = np.ones((seq_count,nr_features))
    # Store contig_lengths of the sequences
    contig_lengths = np.zeros(seq_count)
    
    contigs_id = []
    for i,seq in enumerate(SeqIO.parse(comp_file,"fasta")):
        contigs_id.append(seq.id)
        contig_lengths[i] = len(seq)
        for kmer_tuple in window(seq.seq.tostring().upper(),kmer_len):
            kmer = "".join(kmer_tuple)
            if kmer in feature_mapping:
                composition[i,feature_mapping[kmer]] += 1

    contig_lengths = p.Series(contig_lengths,index=contigs_id,dtype=float)
    composition = p.DataFrame(composition,index=contigs_id,dtype=float)

    # Select contigs to cluster on, namely the sequences longer than the threshold.
    # The total kmer count without pseudo counts is related to the sequence length through:
    # 
    # Kmer_count = Seq_length - kmer_len + 1
    threshold_filter = contig_lengths > threshold
    
    #log(p_ij) = log[(X_ij +1) / rowSum(X_ij+1)]
    composition = np.log(composition.divide(composition.sum(axis=1),axis=0))
    
    logging.info('Successfully loaded composition data.')
    return composition,contig_lengths,threshold_filter

def load_coverage(cov_file, contig_lengths, no_cov_normalization, add_total_coverage=False, read_length=100):
    #Coverage import, file has header and contig ids as index
    cov = p.read_table(cov_file, header=0, index_col=0)

    temp_cov_range = None
    # cov_range variable left here for historical reasons. Can be removed entirely
    cov_range = (cov.columns[0],cov.columns[-1])

    # Adding pseudo count
    cov.ix[:,cov_range[0]:cov_range[1]] = cov.ix[:,cov_range[0]:cov_range[1]].add(
            (read_length/contig_lengths),
            axis='index')

    if not no_cov_normalization:
        #Normalize per sample first
        cov.ix[:,cov_range[0]:cov_range[1]] = \
            _normalize_per_sample(cov.ix[:,cov_range[0]:cov_range[1]])

    # Total coverage should be calculated after per sample normalization
    if add_total_coverage:
        cov['total_coverage'] = cov.ix[:,cov_range[0]:cov_range[1]].sum(axis=1)
        temp_cov_range = (cov_range[0],'total_coverage')
    
    if not no_cov_normalization:
        # Normalize contigs next
        cov.ix[:,cov_range[0]:cov_range[1]] = \
            _normalize_per_contig(cov.ix[:,cov_range[0]:cov_range[1]])

    if temp_cov_range:
        cov_range = temp_cov_range

    # Log transform
    cov.ix[:,cov_range[0]:cov_range[1]] = np.log(
        cov.ix[:,cov_range[0]:cov_range[1]])

    logging.info('Successfully loaded coverage data.')
    return cov,cov_range
    
def _normalize_per_sample(arr):
    """ Divides respective column of arr with its sum. """
    return arr.divide(arr.sum(axis=0),axis=1)

def _normalize_per_contig(arr):
    """ Divides respective row of arr with its sum. """
    return arr.divide(arr.sum(axis=1),axis=0)
    

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        kmer = ''.join(kmer)
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = ''.join([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in xrange(i):
            next(el, None)
    return izip(*els)
