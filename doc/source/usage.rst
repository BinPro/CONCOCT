
Basic Usage
===========

Cut contigs into smaller parts::

    cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa


Generate table with coverage depth information per sample and subcontig.
This step assumes the directory 'mapping' contains sorted and indexed bam files where each sample has been mapped against the original contigs::

    concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv


Run concoct::

    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/


Merge subcontig clustering into original contig clustering::

    merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv


Extract bins as individual FASTA::

    mkdir concoct_output/fasta_bins
    extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
