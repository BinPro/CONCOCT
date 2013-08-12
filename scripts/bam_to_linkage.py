#!/usr/bin/env python
"""
@author: inodb
"""
import os
import argparse

from Bio import SeqIO
from Bio.SeqUtils import GC
import pysam


#TODO: same function as in gen_input_table.py
def get_gc_and_len_dict(fastafile):
    """Creates a dictionary with the fasta id as key and GC and length as keys
    for the inner dictionary."""
    out_dict = {}

    for rec in SeqIO.parse(fastafile, "fasta"):
        out_dict[rec.id] = {}
        out_dict[rec.id]["GC"] = GC(rec.seq)
        out_dict[rec.id]["length"] = len(rec.seq)

    return out_dict


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


def get_orientation(read, regionlength):
    """Determine orientation of pair."""
    if pair_is_inward(read, regionlength):
        return "inward"
    elif pair_is_outward(read, regionlength):
        return "outward"
    else:
        return "inline"


def get_linkage_info_dict_at_loc(bamh, chromosome, start, end, regionlength, samplename, out_dict={}):
    """Returns a linkage information dictionary given a specific location on a
    contig. The linkage information dictionary shows which contigs are linked
    to what other contigs by pairs. The pairs can be inward pointing, outward
    or in line. The outer dictionary contains contig names, as well as the
    inner to indicate linkage between the two contigs. The information is
    stored twice (for both contigs). The 2nd inner dictionary has a sample name
    and the third is a counter for inner/outward/inline pairs that link the two
    contigs."""
    for read in bamh.fetch(chromosome, start, end):
        # check if pair links to another contig
        # only look at pair1, to prevent counting links twice
        if read.is_paired and read.tid != read.mrnm and read.is_read1:
            ref = bamh.getrname(read.tid)
            refm = bamh.getrname(read.mrnm)

            # Init with empty dic if contig is not in there yet
            if ref not in out_dict:
                out_dict[ref] = {refm : {}}
            elif refm not in out_dict[ref]:
                out_dict[ref][refm] = {}
            # Also check if refm is in the dict, because it could already have
            # links to another contig.
            if refm not in out_dict:
                out_dict[refm] = {ref : {}}
            elif ref not in out_dict[refm]:
                out_dict[refm][ref] = {}

            # Add one link
            ori = get_orientation(read, regionlength)
            try:
                d = out_dict[ref][refm][samplename]
            except KeyError:
                # Sample hasn't been added yet
                out_dict[ref][refm][samplename] = {"inward" : 0, "outward" : 0, "inline" : 0}
                out_dict[refm][ref][samplename] = {"inward" : 0, "outward" : 0, "inline" : 0}

            out_dict[ref][refm][samplename][ori] += 1
            out_dict[refm][ref][samplename][ori] += 1

    return out_dict


def get_linkage_info_dict_one_fetch(bamfile, regionlength, samplename, out_dict={}):
    """Creates a two-dimensional dictionary of linkage information between
    contigs. Fetch is called only once on the bam file, which is faster than a
    fetch for each contig individually. About twice as fast for 100 reads."""
    bamh = pysam.Samfile(bamfile)

    for read in bamh.fetch():
        # check if pair links to another contig
        # only look at pair1, to prevent counting links twice
        # check if both contigs are bigger than twice the regionlength
        # (arbitrary contig length cut-off)
        # check if pairs are both at one end of the contig
        if read.is_paired \
          and read.tid != read.mrnm \
          and read.is_read1 \
          and bamh.lengths[read.tid] > 2 * regionlength \
          and bamh.lengths[read.mrnm] > 2 * regionlength \
          and (read.pos <= regionlength or read.pos >= bamh.lengths[read.tid] - regionlength) \
          and (read.mpos <= regionlength or read.mpos >= bamh.lengths[read.mrnm] - regionlength):
            ref = bamh.getrname(read.tid)
            refm = bamh.getrname(read.mrnm)

            # Init with empty dic if contig is not in there yet
            if ref not in out_dict:
                out_dict[ref] = {refm : {}}
            elif refm not in out_dict[ref]:
                out_dict[ref][refm] = {}
            # Also check if refm is in the dict, because it could already have
            # links to another contig.
            if refm not in out_dict:
                out_dict[refm] = {ref : {}}
            elif ref not in out_dict[refm]:
                out_dict[refm][ref] = {}

            # Add one link
            ori = get_orientation(read, regionlength)
            try:
                d = out_dict[ref][refm][samplename]
            except KeyError:
                # Sample hasn't been added yet
                out_dict[ref][refm][samplename] = {"inward" : 0, "outward" : 0, "inline" : 0}
                out_dict[refm][ref][samplename] = {"inward" : 0, "outward" : 0, "inline" : 0}

            out_dict[ref][refm][samplename][ori] += 1
            out_dict[refm][ref][samplename][ori] += 1

    return out_dict


def get_linkage_info_dict(fdict, bamfile, regionlength, samplename, out_dict={}):
    """Creates a two-dimensional dictionary of linkage information between
    contigs."""
    bamh = pysam.Samfile(bamfile)

    for contig in fdict:
        # Check ends if 2 * regionlength smaller than the contig length
        # skip shorter contigs for linkage
        if regionlength * 2 < fdict[contig]["length"]:
            out_dict = get_linkage_info_dict_at_loc(bamh, contig, 0, regionlength, regionlength,
                samplename, out_dict)
            out_dict = get_linkage_info_dict_at_loc(bamh, contig,
                fdict[contig]["length"] - regionlength - 1, fdict[contig]["length"], regionlength,
                samplename, out_dict)

    return out_dict


def print_linkage_info(fastafile, bamfiles, samplenames, regionlength):
    """Prints a linkage information table. Format as

    contig1<TAB>contig2<TAB>nr_links_inward_n<TAB>nr_links_outward_n

    where n represents sample name. Number of columns is thus 2 + 2 * n.
    """
    # With one fetch determining length beforehand is no longer necessary
    #fdict = get_gc_and_len_dict(fastafile)

    assert(len(bamfiles) == len(samplenames))

    linkdict = {}
    for i in range(len(bamfiles)):
        linkdict = get_linkage_info_dict_one_fetch(bamfiles[i], regionlength,
            samplenames[i], linkdict)

    # Header
    print ("%s\t%s" + "\t%s" * len(bamfiles)) % (("contig1", "contig2") +
        tuple(["nr_links_inward_%s\tnr_links_outward_%s" % (s, s) for s in samplenames]))

    # Content
    for contig in linkdict:
        for contig2 in linkdict[contig]:
            nrlinksl = []

            for s in samplenames:
                if s in linkdict[contig][contig2]:
                    nrlinksl = nrlinksl + [linkdict[contig][contig2][s]["inward"], linkdict[contig][contig2][s]["outward"]]
                else:
                    nrlinksl = nrlinksl + [0, 0]

            if sum(nrlinksl) > 0:
                print ("%s\t%s" + "\t%i" * 2 * len(bamfiles)) % ((contig, contig2) + tuple(nrlinksl))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="Contigs fasta file")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr as bamfiles.")
    parser.add_argument("--regionlength", type=int, default=500, help="Linkage is checked "
        "on both ends of the contig, this parameter specifies the search region "
        "length. (500)")
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

    print_linkage_info(args.fastafile, args.bamfiles, samplenames, args.regionlength)
