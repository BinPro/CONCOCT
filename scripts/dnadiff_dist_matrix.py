#!/usr/bin/env python
"""
Output distance matrix between fasta files using dnadiff from MUMmer.
"""
import argparse
import subprocess
import utils
import re
from os.path import join as ospj
import sys
import numpy as np


class CmdException(Exception):
    """Exception for a shell command that did not exit with returncode 0."""
    def __init__(self, cmd, stdout, stderr, returncode):
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode

    def __str__(self):
        return \
          "cmd:\n{}\n" \
          "returncode:\n{}\n" \
          "stdout:\n{}\n" \
          "stderr:\n{}\n".format(self.cmd, self.returncode, self.stdout,
          self.stderr)


class MUMmerReport(object):
    """Represents .report file from MUMmer's dnadiff. Stores TotalBases and
    AlignedBases (%) stats."""
    def __init__(self, report_file):
        self.report_file = report_file

        with open(report_file) as f:
            one_to_one_parsed = False

            for line in f:
                if line.startswith("TotalBases"):
                    self.tot_bases = [int(b) for b in line.split()[1:]]
                if line.startswith("AlignedBases"):
                    # store percentage
                    self.aligned_bases = [float(p) for p in
                            re.findall(r'\((.*?)\%\)', line)]
                if not one_to_one_parsed and line.startswith("AvgIdentity"):
                    self.avg_identity = [float(p) for p in line.split()[1:]]
                    one_to_one_parsed = True


def run_dnadiff(fasta1, fasta2, prefix):
    """Runs MUMmer's dnadiff"""
    cmd = "dnadiff -p {prefix} {f1} {f2}".format(f1=fasta1, f2=fasta2,
            prefix=prefix)
    with open("{}.cmd".format(prefix), "w") as f:
        f.write(cmd)
    proc = subprocess.Popen(cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True)
    stdout, stderr = proc.communicate()
    if proc.returncode > 0:
        raise CmdException(cmd, stdout, stderr, proc.returncode)


def run_dnadiff_pairwise(fasta_files, fasta_names, output_folder):
    """Runs MUMmer's dnadiff pairwise for given fasta_files. Uses fasta_names
    to organize output folders for dnadiff as fastaname1_vs_fastaname2."""
    assert len(fasta_files) == len(fasta_names)

    for i in range(len(fasta_files)):
        for j in range(i + 1, len(fasta_files)):
            out_dir = ospj(output_folder, "{fn1}_vs_{fn2}".format(
                fn1=fasta_names[i], fn2=fasta_names[j]))
            utils.mkdir_p(out_dir)
            run_dnadiff(fasta_files[i], fasta_files[j], ospj(out_dir, "out"))


def get_dist_matrix(pairwise_folder, fasta_names, min_coverage):
    """Returns distance matrix from folder constructed with
    run_dnadiff_pairwise"""
    matrix = np.array(len(fasta_names) * [len(fasta_names) * [0.0]])
    for i in range(len(fasta_names)):
        for j in range(i + 1, len(fasta_names)):
            repfile = ospj(pairwise_folder, "{fn1}_vs_{fn2}".format(
                fn1=fasta_names[i], fn2=fasta_names[j]), "out.report")
            mumr = MUMmerReport(repfile)
            if min(mumr.aligned_bases) >= min_coverage:
                # take distance as 1 - AvgIdentity (same for both ref and qry)
                matrix[i][j] = 1.0 - (mumr.avg_identity[0] / 100.0)
            else:
                matrix[i][j] = 1.0
            matrix[j][i] = matrix[i][j]

    return matrix


def plot_dist_matrix(matrix, fasta_names, output_file):
    """Cluster the distance matrix hierarchically and plot using seaborn.
    Average linkage method is used."""
    import seaborn as sns
    import pandas as pd

    pdm = pd.DataFrame(matrix, index=fasta_names, columns=fasta_names)
    clustergrid = sns.clustermap(pdm, method='average')
    clustergrid.savefig(output_file)


def parse_input():
    """Return input arguments using argparse"""
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_folder", help="Output folder")
    parser.add_argument("fasta_files", nargs='+',
            help="fasta files to compare pairwise using MUMmer's dnadiff")
    parser.add_argument("--min_coverage", type=float, default=50,
            help="Minimum coverage of bin in percentage to calculate distance "
                 "otherwise distance is 1.")
    parser.add_argument("--fasta_names", default=None, help="File with names "
            "for fasta file, one line each. Could be sample names, bin names, "
            "genome names, whatever you want. The names are used when storing "
            "the MUMmer dnadiff results as in "
            "output_folder/fastaname1_vs_fastaname2/")
    parser.add_argument("--hclust_plot_file", default=None, help="Cluster the "
            "distance matrix hierarchically and plot using seaborn to given "
            "file. Uses average linkage method for clustering. Type of image is "
            "determined by extension.")
    parser.add_argument("--skip_dnadiff", action="store_true", help="Skips "
            "running MUMmer and uses output_folder as given input to "
            "calculate the distance matrix. Expects dnadiff output as "
            "output_folder/fastaname1_vs_fastaname2/out.report")
    parser.add_argument("--skip_matrix", action="store_true", help="Skips "
            "Calculating the distance matrix i.e. just runs dnadiff.")
    args = parser.parse_args()
    # Get fasta names
    if args.fasta_names is not None:
        fasta_names = [s[:-1] for s in open(args.fasta_names).readlines()]
        if len(fasta_names) != len(args.fasta_files):
            raise Exception("Nr of names in fasta_names should be equal to nr "
                    "of given fasta_files")
    else:
        fasta_names = list(range(len(args.fasta_files)))
    if args.skip_dnadiff and args.skip_matrix:
        raise Exception("If running dnadiff and calculating the distance matrix "
                "are both skipped, the program does not run any steps at all.")

    return args.output_folder, args.fasta_files, fasta_names, \
           args.min_coverage, args.skip_dnadiff, args.skip_matrix, \
           args.hclust_plot_file


def main(output_folder, fasta_files, fasta_names, min_coverage, skip_dnadiff=False, skip_matrix=False, hclust_plot_file=None):
    """Output distance matrix between fasta files using MUMmer's dnadiff"""
    if not skip_dnadiff:
        run_dnadiff_pairwise(fasta_files, fasta_names, output_folder)
    if not skip_matrix:
        matrix = get_dist_matrix(output_folder, fasta_names, min_coverage)
        np.savetxt(sys.stdout, matrix, fmt="%.2f", delimiter="\t")
    if hclust_plot_file:
        plot_dist_matrix


if __name__ == "__main__":
    main(*parse_input())
