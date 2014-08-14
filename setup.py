#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
from distutils.core import Extension
from setuptools.command.install import install as _install

version = '0.3.2'

def _post_install():
    from subprocess import call
    call(['echo','####### Will now install bash completion for python scripts'])
    call(['activate-global-python-argcomplete', '--user'])
    call(['echo','####### To activate bash completion for Concoct, you need ',
          'to add the following line in your ~/.bashrc file: \n'])
    call(['echo','source ~/.bash_completion.d/python-argcomplete.sh\n'])
    call(['echo','####### If this gives you errors in bash, you can try ',
          'the following instead'])
    call(['echo','eval "$(register-python-argcomplete concoct)"\n'])

class install(_install):
    def run(self):
        _install.run(self)
        self.execute(_post_install,(),
             msg="\n############ Running post install task ##############\n")

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
      cmdclass={'install':install},
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
                        'argcomplete>=0.6.7',
                        'nose==1.3.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
