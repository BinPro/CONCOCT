#!/usr/bin/env python
"""
@author: inodb
"""
import os
import argparse
import multiprocessing

import pysam

import ipdb


class Read:
    def __init__(self):
        self.align_count = 0

    def count(self):
        self.align_count += 1

    def get_value(self):
        return 1 / self.align_count


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
            return "inline"
        else:
            return "inward_or_outward"
    else:
        return get_orientation_tips(read, regionlength)


def get_orientation_tips(read, regionlength):
    """Determine orientation of pair if the reads are located on the tips of
    the contig."""
    if pair_is_inward(read, regionlength):
        return "inward"
    elif pair_is_outward(read, regionlength):
        return "outward"
    else:
        return "inline"


def is_within_region(readpos, contiglength, readlength, regionlength):
    return readpos - readlength <= regionlength \
      or readpos >= contiglength - regionlength


def is_link(read):
    return read.is_paired \
      and read.tid != read.mrnm


def get_linkage_info_dict_one_fetch(bamfile, readlength, min_contig_length, regionlength, fullsearch):
    """Creates a two-dimensional dictionary of linkage information between
    contigs. Fetch is called only once on the bam file."""
    bamh = pysam.Samfile(bamfile)
    linkdict = {}
    read_count_dict = {}

    for read in bamh:
        # check if contig meets length threshold
        # check if the read falls within specified region
        if bamh.lengths[read.tid] >= min_contig_length \
          and (fullsearch
          or is_within_region(read.pos, bamh.lengths[read.tid], readlength, regionlength)):
            ref = bamh.getrname(read.tid)
            read_count_dict[ref] = read_count_dict.get(ref, 0) + 1

            # look at first read to prevent counting links twice
            # linked contig should be within threshold as well
            # linked mate should be within specified region
            if read.is_read1 \
              and read.is_paired \
              and bamh.lengths[read.mrnm] >= min_contig_length \
              and (fullsearch
              or is_within_region(read.mpos, bamh.lengths[read.mrnm], readlength, regionlength)):
                refm = bamh.getrname(read.mrnm)

                #if refm == "contig-2482000051":
                #    ipdb.set_trace()

                init_link_row = lambda: {"inward" : 0, "outward" : 0, "inline" : 0, "inward_or_outward" : 0}
                # Init if contig is not in there yet
                if ref not in linkdict:
                    linkdict[ref] = {refm : init_link_row()}
                elif refm not in linkdict[ref]:
                    linkdict[ref][refm] = init_link_row()
                # Also check if refm is in the dict, because it could already have
                # links to another contig.
                if refm not in linkdict:
                    linkdict[refm] = {ref : init_link_row()}
                elif ref not in linkdict[refm]:
                    linkdict[refm][ref] = init_link_row()

                # Add one link
                ori = get_orientation(read, regionlength, readlength, bamh.lengths[read.tid], bamh.lengths[read.mrnm])
                linkdict[ref][refm][ori] += 1
                linkdict[refm][ref][ori] += 1

    return linkdict, read_count_dict


def parallel_get_linkage_info_dict(args):
    i, bamfile, readlength, min_contig_length, regionlength = args
    return i, get_linkage_info_dict_one_fetch(bamfile, readlength, min_contig_length, regionlength)


def print_linkage_info(fastafile, bamfiles, samplenames, readlength, min_contig_length, regionlength, max_n_processors, fullsearch):
    """Prints a linkage information table. Format as

    contig1<TAB>contig2<TAB>nr_links_inward_n<TAB>nr_links_outward_n

    where n represents sample name. Number of columns is thus 2 + 2 * n.
    """
    assert(len(bamfiles) == len(samplenames))

    # Determine links in parallel
    linkdict_args = []
    for i in range(len(bamfiles)):
        linkdict_args.append((i, bamfiles[i], readlength, min_contig_length, regionlength))

    #n_processes = min(multiprocessing.cpu_count(), max_n_processors)
    #pool = multiprocessing.Pool(processes=n_processes)
    #poolrv = pool.map(parallel_get_linkage_info_dict, linkdict_args)

    ## order parallel link results
    #linkdict = {}
    #for rv in poolrv:
    #    linkdict[samplenames[rv[0]]] = rv[1][0]
    linkdict = {}
    read_count_dict = {}
    for s in samplenames:
        linkdict[s], read_count_dict[s] = get_linkage_info_dict_one_fetch(bamfiles[0], readlength, min_contig_length, regionlength, fullsearch)

    # Header
    print ("%s\t%s" + "\t%s" * len(bamfiles)) % (("contig1", "contig2") +
        tuple(["nr_links_inward_%s\tnr_links_outward_%s\tnr_links_inline_%s\tnr_links_inward_or_outward_%s" % (s, s, s, s) for s in samplenames]))

    # Content
    allcontigs = set([k for k in linkdict[s].keys() for s in samplenames])
    for c in allcontigs:
        for c2 in allcontigs:
            nrlinksl = []

            for s in samplenames:
                if c in linkdict[s] and c2 in linkdict[s][c]:
                    nrlinksl = nrlinksl + [linkdict[s][c][c2]["inward"],
                                           linkdict[s][c][c2]["outward"],
                                           linkdict[s][c][c2]["inline"],
                                           linkdict[s][c][c2]["inward_or_outward"]]
                else:
                    nrlinksl = nrlinksl + [0, 0, 0, 0]

            if sum(nrlinksl) > 0:
                print ("%s\t%s" + "\t%i" * 4 * len(bamfiles)) % ((c, c2) + tuple(nrlinksl))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="Contigs fasta file")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr as bamfiles.")
    parser.add_argument("--regionlength", type=int, default=500, help="Linkage "
            "is checked on both ends of the contig. This parameter specifies the "
            "search region length in bases, e.g. setting this to 500 means "
            "check for linkage within 500 bases on both ends. (500)")
    parser.add_argument("--fullsearch", action='store_true')
    parser.add_argument('-m', '--max_n_processors', type=int, default=1,
        help="Specify the maximum number of processors to use, if absent, one processor will be used.")
    parser.add_argument("--readlength", type=int, default=100, help="Specify untrimmed read length of reads.")
    parser.add_argument("--mincontiglength", type=int, default=0, help="Length threshold for considered contigs.")

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

    print_linkage_info(args.fastafile, args.bamfiles, samplenames, args.readlength, args.mincontiglength, args.regionlength, args.max_n_processors, args.fullsearch)
