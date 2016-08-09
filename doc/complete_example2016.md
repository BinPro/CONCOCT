#Complete Example V0.3 STAMPS 2016
=================================
This documentation page aims to be a complete example walk through for the usage of the CONCOCT package version 0.3 except for the assembly and COG annotations step.

It is not required to run all steps. The output files for each step are in the test data repository. At the end of this example the results should be the same as the results in the corresponding test data repository: https://github.com/BinPro/CONCOCT-test-data/releases. The version numbers listed above are the ones used to generate the results in that repository. Using newer versions will probably not be a problem, but your results may be different in that case.

##Login to the class servers
-----------------------

    ssh yourname@class.mbl.edu
    ssh yourname@classxx


##Setting up the test environment
-------------------------------
Move to your home directory, create a folder where you want all the output from this example to go:

```
    cd ~
    mkdir CONCOCT-complete-example
    cd CONCOCT-complete-example
```

Set three variables with full paths. One pointing to the root directory of the ```CONCOCT``` software, one pointing to the test data repository, named ```CONCOCT_TEST``` and one to the directory we just created. 

```
    export CONCOCT=/class/stamps-software/CONCOCT
    export CONCOCT_TEST=/class/stamps-software/CONCOCT-test-data
    export CONCOCT_EXAMPLE=$HOME/CONCOCT-complete-example
```

You can see the full path of a directory you are located in by running the command ```pwd```.

Copy in the example data from the class folder:
```
    cp $CONCOCT_TEST/Example.tar.gz .
```

And extract:
```
    tar -xvzf Example.tar.gz
```

##Assembling Metagenomic Reads
----------------------------
The first step in the analysis is to assemble all reads into contigs, even 
for this small example, 1 million reads per sample and 16 samples, this is computationally 
intensive. There are a number of assemblers available if you want a quick results \megahit generally performs well but Spades is also good. The best choice is data set dependent.

Then assemble the reads. We recommend megahit for this. To perform the assembly I ran the following commands, please **do not run this**:

```
cd Example 
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 8 -o Assembly --presets meta > megahit.out&
```

However, I **do not** suggest you do this now instead copy the Assembly directory from the class folders:

```
    cp -r $CONCOCT_TEST/Assembly .
```


##Cutting up contigs
----------------------------
In order to give more weight to larger contigs and mitigate the effect of assembly errors we cut up the contigs into chunks of 10 Kb. The final chunk is appended to the one before it if it is < 10 Kb to prevent generating small contigs. This means that no contig < 20 Kb is cut up. We use the script ``cut_up_fasta.py`` for this:

Now lets cut up the contigs and index for the mapping program bwa:
```
mkdir contigs
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m Assembly/final.contigs.fa > contigs/final_contigs_c10K.fa
```


##Map the Reads onto the Contigs
------------------------------
After assembly we map the reads of each sample back to the assembly using [bwa](https://github.com/lh3/bwa).

We are not going to perform mapping ourselves. Instead just copy the pre-calculated files:

```
    cp -r $CONCOCT_TEST/map .
```

These are the commands we ran (**do not run this**).

```
mkdir Map

for file in Example/*R1.fastq
do 
   
   stub=${file%_R1.fastq}
   stub2=${stub#Example\/}	
   echo $stub

   file2=${stub}_R2.fastq

   bwa mem -t 4 contigs/final_contigs_c10K.fa $file $file2 > Map/${stub2}.sam
done
```

Followed by these (** do not run**):

```
for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub	
    samtools view -h -b -S $file > ${stub}.bam
    samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam
    samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam
    bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g contigs/final_contigs_c10K.len > ${stub}_cov.txt
done
```

##Generate coverage table
------------------------

Given the processed mapping files we can now generate the coverage table:

```
$CONCOCT/scripts/Collate.pl Map | tr "," "\t" > Coverage.tsv
```

This is a simple tab delimited file with the coverage of each contig per sample.


##Run concoct
-----------

To see possible parameter settings with a description run

```
    concoct --help
```

We will only run concoct for some standard settings here. The only one we vary is the cluster number which should be at least twice the number of genomes in your 
co-assembly (see discussion below of how to estimate this). In this case we know it is 
around 20 so run concoct with 40 as the maximum number of cluster `-c 40`:

```
mkdir Concoct
cd Concoct
mv ../Coverage.tsv .
concoct --coverage_file Coverage.tsv --composition_file ../contigs/final_contigs_c10K.fa
cd ..
```

When concoct has finished the message "CONCOCT Finished, the log shows how it went." is piped to stdout. The program generates a number of files in the output directory that can be set with the `-b` parameter and will be the present working directory by default. 

##Contig annotation

We are going to annotate COGs on our contigs. You first need to find genes on the contigs and functionally annotate these. Here we used prodigal (https://github.com/hyattpd/Prodigal) for gene prediction and annotation, but you can use anything you want (**do not run this**):

```
mkdir Annotate_gt1000
cd Annotate_gt1000
python $CONCOCT/scripts/LengthFilter.py -m 1000 ../contigs/final_contigs_c10K.fa > final_contigs_gt1000_c10K.fa
prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff
```

Then we assigned COGs (but on another server) - **do not run this**

```
export COGSDB_DIR=/home/opt/rpsblast_db
nohup $CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 32 -r 1 > r.out&
```

Instead just copy whole directory into your example dir (do run this)

```
cd $CONCOCT_EXAMPLE
cp -r $CONCOCT_TEST/Annotate_gt1000 .
```


##Evaluate output
---------------

This will require that you have Rscript with the R packages [gplots](http://cran.r-project.org/web/packages/gplots/index.html), [reshape](http://cran.r-project.org/web/packages/reshape/index.html), [ggplot2](http://cran.r-project.org/web/packages/ggplot2/index.html), [ellipse](http://cran.r-project.org/web/packages/ellipse/index.html), [getopt](http://cran.r-project.org/web/packages/getopt/index.html) and [grid](http://cran.r-project.org/web/packages/grid/index.html) installed. The package grid does not have to be installed for R version > 1.8.0

First we can visualise the clusters in the first two PCA dimensions:

```
cd $CONCOCT_EXAMPLE
mkdir evaluation-output
Rscript $CONCOCT/scripts/ClusterPlot.R -c Concoct/clustering_gt1000.csv -p Concoct/PCA_transformed_data_gt1000.csv -m Concoct/pca_means_gt1000.csv -r Concoct/pca_variances_gt1000_dim -l -o evaluation-output/ClusterPlot.pdf
```

<figs/ClusterPlot.pdf>

To visualise your plots you will have to copy them off the server to a local directory. Move to that local directory on your computer and type:

    scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/ClusterPlot.pdf .

We can also compare the clustering to species labels. For this test data set we know these labels, they are given in the file ```$CONCOCT_TEST/AssignGenome/clustering_gt1000_smap.csv```. For real data labels may be obtained through taxonomic classification.
In either case we provide a script Validate.pl for computing basic metrics on the cluster quality. Lets copy the validation data into our directory:

```
cd $CONCOCT_EXAMPLE
cp -r $CONCOCT_TEST/AssignGenome .
```

And run our validation script:
```
cd evaluation-output
$CONCOCT/scripts/Validate.pl --cfile=../Concoct/clustering_gt1000.csv --sfile=clustering_gt1000_smap.csv --ffile=../Annotate_gt1000/final_contigs_gt1000_c10K.fa
```


    cd $CONCOCT_EXAMPLE
    cp $CONCOCT_TEST/evaluation-output/clustering_gt1000_s.csv evaluation-output/
    $CONCOCT/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv --ffile=contigs/velvet_71_c10K.fa
    

This script requires the clustering output by concoct ```concoct-output/clustering_gt1000.csv``` these have a simple format of a comma separated file listing each contig id followed by the cluster index and the species labels that have the same format but with a text label rather than a cluster index. The script should output:

    N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    684	684	6.8023e+06	5	4	0.897224	0.999604	0.841911	0.911563	0.823200


This gives the no. of contigs N clustered, the number with labels M, the number of unique labels S, the number of clusters K, the recall, the precision, the normalised mutual information (NMI), the Rand index, and the adjusted Rand index. It also generates a file called a `confusion matrix` with the frequencies of each species in each cluster. We provide a further script for visualising this as a heatmap:

    Rscript $CONCOCT/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf

This generates a file with normalised frequencies of contigs from each cluster across species:

<https://github.com/BinPro/CONCOCT-test-data/tree/master/evaluation-output/clustering_gt1000_conf.pdf>

To view this also download off server to your local directory:

    scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/clustering_gt1000_conf.pdf .

Validation using single-copy core genes
---------------------------------------

We can also evaluate the clustering based on single-copy core genes. You first need to find genes on the contigs and functionally annotate these. Here we used prodigal (https://github.com/hyattpd/Prodigal) for gene prediction and annotation, but you can use anything you want:

    cd $CONCOCT_EXAMPLE
    mkdir -p $CONCOCT_EXAMPLE/annotations/proteins
    prodigal -a annotations/proteins/velvet_71_c10K.faa \
             -i contigs/velvet_71_c10K.fa \
             -f gff -p meta  > annotations/proteins/velvet_71_c10K.gff

We used RPS-Blast to COG annotate the protein sequences using the script ``RSBLAST.sh``. This is not available on the assemblers so we will just copy the output:
~~You need to set the evironmental variable ``COGSDB_DIR``:~~

~~export COGSDB_DIR=/proj/b2010008/nobackup/database/cog_le/~~
    
~~The script furthermore requires GNU parallel and rpsblast. Here we run it on eight cores:~~

~~$CONCOCT/scripts/RPSBLAST.sh -f annotations/proteins/velvet_71_c10K.faa -p -c 8 -r 1~~
~~mkdir $CONCOCT_EXAMPLE/annotations/cog-annotations~~
~~mv velvet_71_c10K.out annotations/cog-annotations/~~

The blast output has been placed in:

    $CONCOCT_TEST/annotations/cog-annotations/velvet_71_c10K.out
    
Copy to your local directory:

    mkdir $CONCOCT_EXAMPLE/annotations/cog-annotations
    cp $CONCOCT_TEST/annotations/cog-annotations/velvet_71_c10K.out $CONCOCT_EXAMPLE/annotations/cog-annotations/velvet_71_c10K.out

Finally, we filtered for COGs representing a majority of the subject to ensure fragmented genes are not over-counted and generated a table of counts of single-copy core genes in each cluster generated by CONCOCT. Remember to use a real email adress, this is supplied since information is fetched from ncbi using their service eutils, and the email is required to let them know who you are.

    cd $CONCOCT_EXAMPLE
    $CONCOCT/scripts/COG_table.py -b annotations/cog-annotations/velvet_71_c10K.out \
    -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt \
    -c concoct-output/clustering_gt1000.csv \
    --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > evaluation-output/clustering_gt1000_scg.tab

The script requires the clustering output by concoct ```concoct-output/clustering_gt1000.csv```, a file listing a set of SCGs (e.g. a set of COG ids) to use ```scgs/scg_cogs_min0.97_max1.03_unique_genera.txt``` and a mapping of Conserved Domain Database ids (https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) to COG ids ``$CONCOCT/scgs/cdd_to_cog.tsv``.
If these protein sequences were generated by Prokka, the names of the contig ids needed to be recovered from the gff file. Since prodigal has been used, the contig ids instead are recovered from the protein ids using a separator character, in which case only the string before (the last instance of) the separator will be used as contig id in the annotation file. In the case of prodigal the separator that should be used is _ and this is the default value, but other characters can be given through the '--separator' argument.

The output file is a tab-separated file with basic information about the clusters (cluster id, ids of contigs in cluster and number of contigs in cluster) in the first three columns, and counts of the different SCGs in the following columns.

This can also be visualised graphically using the R script:

    cd $CONCOCT_EXAMPLE
    Rscript $CONCOCT/scripts/COGPlot.R -s evaluation-output/clustering_gt1000_scg.tab -o evaluation-output/clustering_gt1000_scg.pdf

To view this also download off server to your local directory:

    scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/clustering_gt1000_scg.pdf .

The plot is downloadable here:

<https://github.com/BinPro/CONCOCT-test-data/tree/master/evaluation-output/clustering_gt1000_scg.pdf>

Incorporating linkage information
---------------------------------

To perform a hierarchical clustering of the clusters based on linkage we simply run:

    $CONCOCT/scripts/ClusterLinkNOverlap.pl --cfile=concoct-output/clustering_gt1000.csv --lfile=concoct-input/concoct_linkage.tsv --covfile=concoct-input/concoct_inputtableR.tsv --ofile=concoct-output/clustering_gt1000_l.csv

The output indicates that the clusters have been reduced from four to three. The new clustering is given by ```concoct-output/clustering_gt1000_l.csv```. This is a significant improvement in recall:

    $CONCOCT/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000_l.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv
    N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    684	684	6.8400e+02	5	3	1.000000	0.997076	0.995805	0.999979	0.999957

The algorithm is explained in more depth in the paper on [arXiv](http://arxiv.org/abs/1312.4038)
