.. CONCOCT documentation master file, created by
   sphinx-quickstart on Wed Jul 23 15:17:29 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CONCOCT's documentation
=======================
CONCOCT "bins" metagenomic contigs. Metagenomic binning is the process of clustering sequences into clusters corresponding to operational taxonomic units of some level.

For any known issues with CONCOCT check the issue tracker:
https://github.com/BinPro/CONCOCT/issues

Features
--------
CONCOCT does unsupervised binning of metagenomic contigs by using nucleotide composition - kmer frequencies - and coverage data for multiple samples. 
CONCOCT can accurately (up to species level) bin metagenomic contigs. For optmial performance:

- Map several samples against your assembled contigs.
- Cut longer contigs into 10 - 20 kb pieces prior to mapping.
- Evaluate your bins using single copy genes.

Installation
------------
For a comprehensive guide on how to install CONCOCT and all its dependencies, see :doc:`installation`.

Contribute
----------

- Issue Tracker: `github <https://github.com/BinPro/CONCOCT/issues>`__
- Source Code: `github <https://github.com/BinPro/CONCOCT>`__

Support
-------
If you are having issues, please let us know.
We have a mailing list located at: concoct-support@lists.sourceforge.net which you can subscribe to `here <https://lists.sourceforge.net/lists/listinfo/concoct-support>`__.

Licence
-------
FreeBSD

Contents:
---------
.. toctree::
   :maxdepth: 2

   self
   installation
   usage
   complete_example
   scripts/index
