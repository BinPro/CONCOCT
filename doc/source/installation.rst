Installation
============

Dependencies
------------

Fundamental dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

::

    python v2.7.*
    gcc
    gsl

These items are prerequisities for the installation of concoct as
described below. The installation procedure varies on different systems,
and described in this README is only how to proceed with a linux
(ubuntu) distribution.

The first item, ``python v2.7.*``, should be installed on a modern
Ubuntu distribution. A c-compiler, e.g. ``gcc``, is needed to compile
the c parts of concoct that uses the GNU Scientific Library ``gsl``. For
linux (ubuntu) this is installed through:

::

    apt-get install build-essential libgsl0-dev

Python packages
~~~~~~~~~~~~~~~

::

    cython>=0.19.2
    numpy>=1.7.1
    scipy>=0.12.0
    pandas>=0.11.0
    biopython>=1.62b
    scikit-learn>=0.13.1

These are the python packages that need to be installed in order to run
concoct. If you follow the installation instructions below, these will
be installed automatically, but are listed here for transparency.

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

-  For assembly, use your favorite, here is one

   -  `Velvet <http://www.ebi.ac.uk/~zerbino/velvet/>`__

      -  In velvet installation directory Makefile, set
         'MAXKMERLENGTH=128', if this value is smaller in the default
         installation.

-  To create the input table (containing average coverage per sample and
   contig)

   -  `BEDTools <https://github.com/arq5x/bedtools2/releases>`__ version
      >= 2.15.0 (only genomeCoverageBed)
   -  `Picard <https://launchpad.net/ubuntu/+source/picard-tools/>`__
      tools version >= 1.110
   -  `samtools <http://samtools.sourceforge.net/>`__ version >= 0.1.18
   -  `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__
      version >= 2.1.0
   -  `GNU parallel <http://www.gnu.org/software/parallel/>`__ version
      >= 20130422
   -  Python packages: ``pysam>=0.6``

-  For validation of clustering using single-copy core genes

   -  `Prodigal <http://prodigal.ornl.gov/>`__ >= 2.60
   -  Python packages: ``bcbio-gff>=0.4``
   -  R packages: ``gplots, reshape, ggplot2, ellipse, getopt`` and
      ``grid``
   -  `BLAST <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/>`__ >= 2.2.28+

If you want to install these dependencies on your own server, you can
take a look at `doc/Dockerfile.all\_dep <doc/Dockerfile.all_dep>`__ for
ideas on how to install them.

Installation
------------

Here we describe two recommended ways of getting concoct to run on your
computer/server. The first option, using Anaconda, should work for any
\*nix (e.g. Mac OS X or Linux) system even where you do not have 'sudo'
rights (e.g. on a common computer cluster). The second option is
suitable for a linux computer where you have root privileges and you
prefer to use a virtual machine where all dependencies to run concoct
are included. Docker does also run on Mac OS X through a virtual machine.
For more information check out the `Docker documentation <http://docs.docker.com/installation/>`__.

Using Anaconda
~~~~~~~~~~~~~~

This instruction shows how to install all dependencies (except the
'Fundamental dependencies' and the 'Optional dependencies' listed above)
using an Anaconda environment. Anaconda is a tool to isolate your python
installation, which allows you to have multiple parallel installations
using different versions of different packages, and gives you a very
convenient and fast way to install the most common scientific python
packages. Anaconda is free but not open source, you can download
Anaconda `here <https://store.continuum.io/cshop/anaconda/>`__.
Installation instructions can be found
`here <http://docs.continuum.io/anaconda/install.html>`__.

After installing Anaconda, create a new environment that will contain
the concoct installation:

::

    conda create -n concoct_env python=2.7.6

After choosing to proceed, run the suggested command:

::

    source activate concoct_env

then install the concoct dependencies into this environment:

::

    conda install cython numpy scipy biopython pandas pip scikit-learn

Finally, download the CONCOCT distribution from
https://github.com/BinPro/CONCOCT/releases (stable) and extract the
files, or clone the repository with github (potentially unstable).
Resolve all dependencies, see above and then execute within the CONCOCT
directory:

::

    python setup.py install

Using Docker
~~~~~~~~~~~~

If you have root access to a machine where you want to install concoct
and you have storage for roughly 2G "virtual machine" then Docker
provides a very nice way to get a Docker image with concoct and its
dependencies installed. This way the only thing you install on your host
system is Docker, the rest is contained in an Docker image. This allows
you to install and run programs in that image without it affecting your
host system. You should `get to know Docker
here <https://docs.docker.com/>`__. You need to `get
Docker installed <https://docs.docker.com/installation/>`__ and
specially if you have
`Ubuntu <http://docs.docker.com/installation/ubuntulinux/>`__.
When Docker is installed you need to download and log into the concoct
image.

We provide a Docker image:

binpro/concoct\_latest contains CONCOCT and all its dependencies for the
:doc:`complete_example` with the exception of
the SCG evaluation.

The following command will then download the image from the Docker image
index, map the Data folder to the image and log you into the docker
image.

::

    sudo docker run -v /home/USER/Data:/opt/Data -i -t binpro/concoct_latest bash

To test concoct you can then do:

::

    $ cd /opt/CONCOCT_latest
    $ nosetests

Which should execute all tests without errors.
