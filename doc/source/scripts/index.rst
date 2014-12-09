CONCOCT Scripts
=================================================
The scripts in the ``CONCOCT/scripts`` directory are not fully maintained. They
implement methods that we apply after binning with CONCOCT. Eventually some of
these methods might make it to a package of their own.

To test all scripts that have tests one could do::

    cd CONCOCT/scripts/tests
    nosetests

Before using a script it would be good to check if its test (in case it has
one) is working for you::

    cd CONCOCT/scripts/tests
    nosetests -s test_script_name

Contents:

.. toctree::
    :maxdepth: 2

    dnadiff_dist_matrix
