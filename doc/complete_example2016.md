#Complete Example V0.3 STAMPS 2016
=================================
This documentation page aims to be a complete example walk through for the usage of the CONCOCT package version 0.3 except for the assembly and COG annotations step.

It is not required to run all steps. The output files for each step are in the test data repository. At the end of this example the results should be the same as the results in the corresponding test data repository. The version numbers listed above are the ones used to generate the results in that repository. Using newer versions will probably not be a problem, but your results may be different in that case.

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


To visualise your plots you will have to copy them off the server to a local directory. Move to that local directory on your computer and type:

    scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/ClusterPlot.pdf .

The figure should look like this:

![Cluster PCA](figs/ClusterPlot.pdf)

We can also compare the clustering to species labels. For this test data set we know these labels, they are given in the file ```$CONCOCT_TEST/AssignGenome/clustering_gt1000_smap.csv```. For real data labels may be obtained through taxonomic classification.
In either case we provide a script Validate.pl for computing basic metrics on the cluster quality. Lets copy the validation data into our directory:

```
cd $CONCOCT_EXAMPLE
cp -r $CONCOCT_TEST/AssignGenome .
```

And run our validation script:
```
cd evaluation-output
$CONCOCT/scripts/Validate.pl --cfile=../Concoct/clustering_gt1000.csv --sfile=../AssignGenome/clustering_gt1000_smap.csv --ffile=../Annotate_gt1000/final_contigs_gt1000_c10K.fa
```

 This script requires the clustering output by concoct ```Concoct/clustering_gt1000.csv``` these have a simple format of a comma separated file listing each contig id followed by the cluster index and the species labels that have the same format but with a text label rather than a cluster index. The script should output:

    N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    9176	9176	5.8291e+07	20	23	0.986723	0.950877	0.976244	0.993650    0.946975


This gives the no. of contigs N clustered, the number with labels M, the number of unique labels S, the number of clusters K, the recall, the precision, the normalised mutual information (NMI), the Rand index, and the adjusted Rand index. It also generates a file called a `confusion matrix` with the frequencies of each species in each cluster. We provide a further script for visualising this as a heatmap:

```
cd $CONCOCT_EXAMPLE
Rscript $CONCOCT/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf
```

This generates a file with normalised frequencies of contigs from each cluster across species:

![Confusion matrix](figs/clustering_gt1000_conf.pdf)


To view this also download off server to your local directory:

```
scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/clustering_gt1000_conf.pdf .
```

##Validation using single-copy core genes
---------------------------------------

We can also evaluate the clustering based on single-copy core genes. We will use our results from above in ```Annotate_gt1000``` to do this.

```
cd $CONCOCT_EXAMPLE 
$CONCOCT/scripts/COG_table.py -b Annotate_gt1000/final_contigs_gt1000_c10K.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c  Concoct/clustering_gt1000.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > evaluation-output/clustering_gt1000_scg.tab
```

The script requires the clustering output by concoct ```concoct-output/clustering_gt1000.csv```, a file listing a set of SCGs (e.g. a set of COG ids) to use ```scgs/scg_cogs_min0.97_max1.03_unique_genera.txt``` and a mapping of Conserved Domain Database ids (https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) to COG ids ``$CONCOCT/scgs/cdd_to_cog.tsv``.

The output file is a tab-separated file with basic information about the clusters (cluster id, ids of contigs in cluster and number of contigs in cluster) in the first three columns, and counts of the different SCGs in the following columns.

This can also be visualised graphically using the R script:

```
cd $CONCOCT_EXAMPLE
Rscript $CONCOCT/scripts/COGPlot.R -s evaluation-output/clustering_gt1000_scg.tab -o evaluation-output/clustering_gt1000_scg.pdf
```

To view this also download off server to your local directory:

```
scp yourname@class.mbl.edu:~/CONCOCT-complete-example/evaluation-output/clustering_gt1000_scg.pdf .
```

The plot is downloadable here:

![scg plot](figs/clustering_gt1000_scg.pdf)

##Comparison to MetaBat

We will also compare to a competitor algorithm released last year [MetaBat](https://bitbucket.org/berkeleylab/metabat) and 
[paper](https://peerj.com/articles/1165/). Note the claim 
"MetaBAT outperforms alternative methods in accuracy and computational efficiency on both synthetic and real metagenome datasets.".

```
module load metabat
mkdir Metabat
cd Metabat
runMetaBat.sh -m 1500 ../contigs/final_contigs_c10K.fa ../Map/Sample*_sub.bam
```

We run it on the same contigs utilising the same bam files. The minimum length that MetaBat will allow is 1500bp. So once this is done to provide a fair comparison we will 
also rerun CONCOCT with this minimum contig length.

```
cd $CONCOCT_EXAMPLE
mkdir Concoct_gt1500
concoct --coverage_file ../Concoct/Coverage.tsv --composition_file ../contigs/final_contigs_c10K.fa -c 40 -l 1500
```

Lets compare the Metabat results with CONCOCT. First we maninpulate the cluster assignments into our format:

```
cd $CONCOCT_EXAMPLE/Metabat
grep ">" final_contigs_c10K.fa.metabat-bins*fa | sed 's/.fa:>/,/g' | sed 's/final_contigs_c10K.fa.metabat-bins-_-m_1500\.//g' | awk -F, '{print $2,$1}' OFS=, > clustering_gt1500.csv
```

Then we run the validation script:
```
 $CONCOCT/scripts/Validate.pl --cfile=clustering_gt1500.csv --sfile=../AssignGenome/clustering_gt1000_smap.csv --ffile=../Annotate_gt1000/final_contigs_gt1000_c10K.fa --ofile=clustering_gt1500_conf.csv
```

    N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
    6649	6649	5.2942e+07	17	29	0.872325	0.999852	0.937318	0.984568	0.868407
    
We can also generate confusion plots and the the SCG table..

```
$CONCOCT/scripts/COG_table.py -b ../Annotate_gt1000/final_contigs_gt1000_c10K.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -cclustering_gt1500.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1500_scg.tab

Rscript $CONCOCT/scripts/COGPlot.R -s clustering_gt1500_scg.tab -o metabat_clustering_gt1500_scg.pdf

Rscript $CONCOCT/scripts/ConfPlot.R  -c clustering_gt1500_conf.csv -o metabat_clustering_gt1500_conf.pdf
```

To generate:

![Metabat scg plot](figs/metabat_clustering_gt1500_scg.pdf)

![Metabat conf plot](figs/metabat_clustering_gt1500_conf.pdf)

