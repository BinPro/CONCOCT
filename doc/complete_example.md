Complete Example
================
This documentation page aims to be a complete example walk through for the usage of the CONCOCT package.


Required software
----------------------
To run the example you need the following software:
* Ray version >= 2.1.0

Downloading test data
-----------------------
Still have to upload the test data somewhere.

Assembling Metagenomic Reads
----------------------------
After obtaining the test data we assemble the contigs with Ray.
Go to the folder with the reads from the test data and run Ray:

    mpiexec -n 1 Ray -k 31 -o out_31 \
        -p Sample118_s1e5_R1.fasta Sample118_s1e5_R2.fasta \
        -p Sample120_s1e5_R1.fasta Sample120_s1e5_R2.fasta \
        -p Sample127_s1e5_R1.fasta Sample127_s1e5_R2.fasta \
        -p Sample134_s1e5_R1.fasta Sample134_s1e5_R2.fasta \
        -p Sample177_s1e5_R1.fasta Sample177_s1e5_R2.fasta \
        -p Sample215_s1e5_R1.fasta Sample215_s1e5_R2.fasta \
        -p Sample230_s1e5_R1.fasta Sample230_s1e5_R2.fasta \
        -p Sample234_s1e5_R1.fasta Sample234_s1e5_R2.fasta \
        -p Sample244_s1e5_R1.fasta Sample244_s1e5_R2.fasta \
        -p Sample261_s1e5_R1.fasta Sample261_s1e5_R2.fasta \
        -p Sample263_s1e5_R1.fasta Sample263_s1e5_R2.fasta \
        -p Sample290_s1e5_R1.fasta Sample290_s1e5_R2.fasta \
        -p Sample302_s1e5_R1.fasta Sample302_s1e5_R2.fasta \
        -p Sample321_s1e5_R1.fasta Sample321_s1e5_R2.fasta \
        -p Sample330_s1e5_R1.fasta Sample330_s1e5_R2.fasta \
        -p Sample343_s1e5_R1.fasta Sample343_s1e5_R2.fasta

Map the Reads onto the Contigs
------------------------------
To create coverage profiles

Generate input files
--------------------
Converg files into the correct format

Run concoct
-----------
To see possible parameter settings with a description run

    concoct --help

We will only run concoct for some standard settings here.
