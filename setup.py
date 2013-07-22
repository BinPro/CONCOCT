#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='concoct',
      version=version,
      description="Clustering contigs with composition and coverage",
      long_description="""\
To be done""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg',
      author_email='brynjar.bjarnason@scilifelab.se',
      url='www.github.com/BinPro/CONCOCT',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["concoct/CONCOCT"],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
