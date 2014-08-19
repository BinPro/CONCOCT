#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
from distutils.core import Extension

version = '0.3.3'

module1 = Extension('vbgmm',
        libraries =['gsl',  'gslcblas'],
        include_dirs = ['c-concoct'],
        sources = ['c-concoct/vbgmmmodule.c'])

setup(name='concoct',
      version=version,
      description="Clustering cONtigs with COverage and ComposiTion",
      long_description="""Concoct is a program that combines three types of 
      information - sequence composition, coverage across multiple sample, 
      and read-pair linkage - to automatically bin metagenomic contigs 
      into genomes. """,
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg, Christopher Quince, Anders Andersson, Ino de Bruijn',
      author_email='binni@binnisb.com',
      maintainer='Johannes Alneberg',
      maintainer_email='johannes.alneberg@scilifelab.se',
      url='https://github.com/BinPro/CONCOCT',
      license='FreeBSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["bin/concoct"],
      include_package_data=True,
      zip_safe=False,
      ext_modules = [module1],
      install_requires=['cython>=0.19.1',
                        'numpy>=1.7.1',
                        'scipy>=0.12.0',
                        'pandas>=0.11.0',
                        'biopython>=1.62b',
                        'scikit-learn>=0.13.1',
                        'nose==1.3.0',
                        'sphinx-rtd-theme>=0.1.6',
                        'Sphinx>=1.2.2'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )

