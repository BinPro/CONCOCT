#!/usr/bin/env python
"""
@author: inodb
"""
import os
import argparse
import multiprocessing

from collections import defaultdict

import pysam


INWARD = 0
OUTWARD = 1
INLINE = 2
INOUTWARD = 3


def read_is_inward(read, regionlength):
    """Determine if read is pointing inward"""
    return (read.pos >= regionlength and not read.is_reverse) or \
           (read.pos <= regionlength and read.is_reverse)


def mate_is_inward(read, regionlength):
    """Determine if mate is pointing inward"""
    return (read.mpos >= regionlength and not read.mate_is_reverse) or \
           (read.mpos <= regionlength and read.mate_is_reverse)


def pair_is_inward(read, regionlength):
    """Determine if pair is pointing inward"""
    return read_is_inward(read, regionlength) and mate_is_inward(read, regionlength)


def pair_is_outward(read, regionlength):
    """Determine if pair is pointing outward"""
    return not read_is_inward(read, regionlength) and not mate_is_inward(read, regionlength)


def is_in_overlapping_region(readpos, contiglength, readlength, regionlength):
    return readpos - readlength <= regionlength and readpos >= contiglength - regionlength


def get_orientation(read, regionlength, readlength, contiglength, contigmlength):
    # Check if one of the reads is in an overlapping region or not within a
    # region. In that case link can only follow two orientations
    # inward_or_outward or inline.
    if is_in_overlapping_region(read.pos, contiglength, readlength, regionlength) \
      or is_in_overlapping_region(read.mpos, contigmlength, readlength, regionlength) \
      or not is_within_region(read.pos, contiglength, readlength, regionlength) \
      or not is_within_region(read.mpos, contigmlength, readlength, regionlength):
        if read.is_reverse == read.mate_is_reverse:
            return INLINE
        else:
            return INOUTWARD
    else:
        return get_orientation_tips(read, regionlength)


def get_orientation_tips(read, regionlength):
    """Determine orientation of pair if the reads are located on the tips of
    the contig."""
    if pair_is_inward(read, regionlength):
        return INWARD
    elif pair_is_outward(read, regionlength):
        return OUTWARD
    else:
        return INLINE


def is_within_region(readpos, contiglength, readlength, regionlength):
    return readpos - readlength <= regionlength \
      or readpos >= contiglength - regionlength


def is_link(read):
    return read.is_paired \
      and read.tid != read.mrnm


def _default_linkdict():
    return defaultdict(_default_link)

def _default_link():
    return [0, 0, 0, 0]

def _default_count():
    return 0

def get_linkage_info_dict(bamfile, readlength, min_contig_length, regionlength, fullsearch):
    """Creates a two-dimensional dictionary of linkage information between
    contigs."""
    linkdict = defaultdict(_default_linkdict)
    read_count_dict = defaultdict(_default_count)

    with pysam.Samfile(bamfile, 'rb') as bamh:
        reflens = bamh.lengths

        for read in bamh:
            # check if read is aligned
            try:
                ref = bamh.getrname(read.tid)
            except ValueError:
                continue

            # check if contig meets length threshold
            # check if the read falls within specified region
            if reflens[read.tid] >= min_contig_length \
              and (fullsearch
              or is_within_region(read.pos, reflens[read.tid], readlength, regionlength)):
                read_count_dict[ref] += 1

                # look at first read to prevent counting links twice
                # linked contig should be within threshold as well
                # linked mate should be within specified region
                if read.is_read1 \
                  and is_link(read) \
                  and reflens[read.mrnm] >= min_contig_length \
                  and (fullsearch
                  or is_within_region(read.mpos, reflens[read.mrnm], readlength, regionlength)):
                    # check if mate is aligned
                    try:
                        refm = bamh.getrname(read.mrnm)
                    except ValueError:
                        continue

                    # Add one link
                    ori = get_orientation(read, regionlength, readlength, reflens[read.tid], reflens[read.mrnm])
                    linkdict[ref][refm][ori] += 1
                    linkdict[refm][ref] = linkdict[ref][refm]

    return dict([(k,v) for k, v in linkdict.iteritems()]), dict([(k, v) for k, v in read_count_dict.iteritems()])


def parallel_get_linkage_info_dict(args):
    i, bamfile, readlength, min_contig_length, regionlength, fullsearch = args
    return i, get_linkage_info_dict(bamfile, readlength, min_contig_length, regionlength, fullsearch)


def print_linkage_info(linkdict, read_count_dict, samplenames):
    """Prints a linkage information table. Format as

    contig1<TAB>contig2<TAB>nr_links_inward_n<TAB>nr_links_outward_n

    where n represents sample name. Number of columns is thus 2 + 2 * n.
    """
    # Header
    print ("%s\t%s" + "\t%s" * len(samplenames)) % (("contig1", "contig2") +
        tuple(["nr_links_inward_%s\tnr_links_outward_%s\tnr_links_inline_%s\tnr_links_inward_or_outward_%s\t"
               "read_count_contig1_%s\tread_count_contig2_%s" % ((s,) * 6) for s in samplenames]))

    # Content
    allcontigs = list(set([k for k in linkdict[s].keys() for s in samplenames]))
    for i in xrange(len(allcontigs)):
        c = allcontigs[i]

        for j in xrange(i + 1, len(allcontigs)):
            c2 =allcontigs[j]

            row = []

            for s in samplenames:
                # check links
                try:
                    row += [
                           linkdict[s][c][c2][INWARD],
                           linkdict[s][c][c2][OUTWARD],
                           linkdict[s][c][c2][INLINE],
                           linkdict[s][c][c2][INOUTWARD],
                           ]
                except KeyError:
                    row += [0, 0, 0, 0]
                # check counts
                try:
                    row += [
                           read_count_dict[s][c],
                           read_count_dict[s][c2]
                           ]
                except KeyError:
                    row += [0, 0]

            # output if contigs have links
            if any([row[i] for i in xrange(6 * len(samplenames)) if i % 6 < 4]):
                print ("%s\t%s" + "\t%i" * 6 * len(samplenames)) % ((c, c2) + tuple(row))


def parse_and_print_linkage_info(fastafile, bamfiles, samplenames, readlength, min_contig_length, regionlength, max_n_cores, fullsearch):
    assert(len(bamfiles) == len(samplenames))

    # Determine links in parallel
    linkdict_args = []
    for i in range(len(bamfiles)):
        linkdict_args.append((i, bamfiles[i], readlength, min_contig_length, regionlength, fullsearch))
    n_processes = min(len(bamfiles), max_n_cores)
    pool = multiprocessing.Pool(processes=n_processes)
    poolrv = pool.map(parallel_get_linkage_info_dict, linkdict_args)

    # order parallel link results
    linkdict = {}
    read_count_dict = {}
    for rv in poolrv:
        linkdict[samplenames[rv[0]]] = rv[1][0]
        read_count_dict[samplenames[rv[0]]] = rv[1][1]
    #for s in samplenames:
    #    linkdict[s], read_count_dict[s] = get_linkage_info_dict(bamfiles[0], readlength, min_contig_length, regionlength, fullsearch)

    print_linkage_info(linkdict, read_count_dict, samplenames)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="Contigs fasta file")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr as bamfiles.")
    parser.add_argument("--regionlength", type=int, default=500, help="Linkage "
            "is checked on both ends of the contig. This parameter specifies the "
            "search region length in bases, e.g. setting this to 500 means "
            "check for linkage within 500 bases on both ends. (500)")
    parser.add_argument("--fullsearch", action='store_true', help="Search "
            "entire contig for links. Orientation is determined differently for "
            "links that fall outside <regionlength>, since one can only "
            "differentiate between pairs being in sequence or not in sequence in "
            "that case.")
    parser.add_argument('-m', '--max_n_cores', type=int, default=multiprocessing.cpu_count(),
        help="Specify the maximum number of cores to use, if absent, all cores "
        "will be used. Each core reads one bamfile, so number of cores used will "
        "never be bigger than the number of bamfiles. (%i)" %
        multiprocessing.cpu_count())
    parser.add_argument("--readlength", type=int, default=100, help="Specify untrimmed read length of reads. (100)")
    parser.add_argument("--mincontiglength", type=int, default=0, help="Length threshold for considered contigs. (0)")

    args = parser.parse_args()

    # Get sample names
    if args.samplenames is not None:
        samplenames = [s[:-1] for s in open(args.samplenames).readlines()]
        if len(samplenames) != len(args.bamfiles):
            raise(Exception("Nr of names in samplenames should be equal to nr "
            "of given bamfiles"))
    else:
        samplenames = [str(i) for i in range(len(args.bamfiles))]

    for bf in args.bamfiles:
        if not os.path.isfile(bf + ".bai"):
            raise(Exception("No index for %s file found, run samtools index "
            "first on bam file." % bf))

    parse_and_print_linkage_info(args.fastafile, args.bamfiles, samplenames, args.readlength, args.mincontiglength, args.regionlength, args.max_n_cores, args.fullsearch)
