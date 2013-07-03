##############################################################################
CALLING CONCOCT
##############################################################################
1.
User calls concoct with:

  contig fasta file
  contig coverage file
  kmer lenght (default 4 and can be omitted from the command line)

same contigs in both files.

2.
For cluster number estimation:
  use default parameters for cluster number (maybe 40 50 60 70 80 ? )
OR
  call concoct with a list of cluster number to try out (if -c is the flag for
  clusters: -c 10 15 20 50 100)
OR
  provide taxonomy file for part of the contigs, count species = X and try
  X*50%, X*75%, X, X*150% X*200%?

3.
Optional how often to run each cluster number with different initialization
default 5

4.
Optional maximum number of iteration per clustering if convergence is not
achieved.
default 100

5.
Output directory where to return output text files. This directory will be
created if it does not exist, and in it there will be created a folder with
the date and time of execution to prevent overriding previous runs.

##############################################################################
THE CLUSTERING
##############################################################################
1.
Read in the fasta file and generate the composition features based on kmer
length.
#!TODO Do we need to do pseudo count or something else with the data?

2.
Generate a boolean array with true if the sum of kmer counts in corresponding
contig > 1000 and false otherwise. We keep the < 1000 kmer count contigs and
will classify them after clustering on > 1000 contigs.

3.
Read the coverage file. It will use the same boolean array as the composition.
Note, we will not assumer the same order, rather both composition and coverage
should use contigs id as an indes to a Pandas DataFrame or Pandas Series. Makes
filtering data based on contig id simple and efficient.
#!TODO Do we need to transform the data?

4.
Run PCA on composition for all the contigs

5.
Select PCA number that explains > 90% of the variance.

6.
Transform both composition and coverage separately to that PCA number

7.
Join the transformed composition and coverage into a single feature vector.
#!TODO if we select PCA number X, the joined vector should be 2*X in size?
#!TODO so we extend composition vector with the coverage vector?

8.
Do scikit-learn Gaussian mixture model with full covariance matrix on the
joined vector for all different cluster sizes and the provided number
of initializations and iterations. We will only cluster the contigs with
True in the boolean vector (contigs > 1000 kmer counts)

##############################################################################
THE RESULTS
##############################################################################
For each cluster number, output the best of the initialized clusters:

1.
Write out the clustering of the > 1000 kmer count to CSV (contigid,clusternumber)

2.
Write out the clusters mean for > 1000 kmer count

3.
Write out the clusters variance for > 1000 kmer count

4.
Write out the clustering of the > 1000 kmer count to CSV (contigid,clusternumber)

5.
Write out the responsibilities for contigs > 1000 kmer count

6.
Write the BIC value for the contigs > 1000 kmer count


7.
Based on current clustering, classify the short contigs to their cluster
 
8.
Write the joined clustering to CSV
