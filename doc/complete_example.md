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
After assembly we map the reads of each sample back to the assembly using bowtie2 and remove PCR duplicates with MarkDuplicates. The coverage histogram for each bam file is computed with bedtools' genomeCoverageBed. The script that calls these programs is provided with CONCOCT. The following command is to be executed in the test data root dir. It creates a folder map and for each sample a subfolder with the alignment against the assembly.

    mkdir map
    cd map
    cp ../out_31/Contigs.fasta raynoscaf_31.fa
    parallel -j 8 mkdir -p {/} '&&' \
        cd {/} '&&' \
        bash $CONCOCT/scripts/map-bowtie2-markduplicates.sh \
            -ct 1 -p '-f' ../{} '$('echo ../{} '|' sed s/R1/R2/')' pair \
            ../raynoscaf_31.fa asm bowtie2 \
        ::: ../*_R1.fasta

The parameters used for `map-bowtie2-markduplicates.sh` are:

* `-c` option to compute coverage histogram with genomeCoverageBed
* `-t` option is number of threads
* `-p` option is the extra parameters given to bowtie2. In this case `-f`.

The five arguments are:
* pair1, the fasta/fastq file with the #1 mates
* pair2, the fasta/fastq file with the #2 mates
* pair_name, a name for the pair used to prefix output files
* assembly, a fasta file of the assembly to map the pairs to
* assembly_name, a name for the assembly, used to postfix outputfiles
* outputfolder, the output files will end up in this folder

Generate coverage table
------------------------
Use the bam files of each sample to create a table with the coverage of each contig per sample. Again executed from the folder `map`:

    python $CONCOCT/scripts/gen_input_table.py --isbedfiles \
        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
        raynoscaf_31.fa */bowtie2/asm_pair-smds.coverage \
    > concoct_inputtable.tsv

Generate linkage table
------------------------
The same bam files can be used to give linkage per sample between contigs:

    python $CONCOCT/scripts/bam_to_linkage.py -m 8 \
        --regionlength 500 --fullsearch \
        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
        raynoscaf_31.fa Sample*/bowtie2/asm_pair-smds.bam \
    > concoct_linkage.tsv

Run concoct
-----------
To see possible parameter settings with a description run

    concoct --help

We will only run concoct for some standard settings here. First we need to parse the input table to just contain the mean coverage for each contig in each sample:

    cut -f1,11-26 concoct_inputtable.tsv > concoct_inputtableR.tsv
