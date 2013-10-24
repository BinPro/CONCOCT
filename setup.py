#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='concoct',
      version=version,
      description="Clustering cONtigs with COverage and ComposiTion",
      long_description="""\
To be done""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg, Christopher Quince, Anders Andersson, Ino de Bruijn',
      author_email='binni@binnisb.com',
      maintainer='Johannes Alneberg',
      maintainer_email='johannes.alneberg@scilifelab.se',
      url='www.github.com/BinPro/CONCOCT',
      license='FreeBSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["concoct/CONCOCT"],
      include_package_data=True,
      zip_safe=False,
      install_requires=['argparse==1.2.1', 'biopython==1.62b',
        'nose==1.3.0', 'numpy==1.7.1', 'pandas==0.11.0',  'scikit-learn==0.13.1',
        'scipy==0.12.0'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
