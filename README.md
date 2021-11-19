# Barcode-counting

The code processes individual timepoint barcode sequencing data from fastq files and assembles them into timecourse trajectories. Data in this form can be used to estimate barcode fitness.

The python code has been broken down in functions that are called at the end.
The functions process the data as follows
Function 'align' merges paired-end reads with pear, then aligns to barcode locus with bowtie2
Function 'barcodes' extracts barcodes and quality scores from the resulting sam files. It stores sequence ID, UMI and sequences and quality scores for BC1 and BC2 in json format
Function 'bartender' collects the BC2 (high complexity) from all libraries in a single file. Then it calls the bartender software in that file and BC2 are clustered
Function 'updateReadsV2' 1. uses known lists of BC1 (Low complexity) to cluster BC1 using levenshtein (or hamming) distance 2. the bartender output to update the BC2 barcodes 3. updates libaray names based on a key
Function 'counts' goes into the updated read files and outputs a file with counts
Function 'aggregateCounts' calculates aggregate counts for assays where the same pools were used
The 'aggregateCounts' output is plotted as histograms and based on that I decide on a filtering step (the minimum number of aggregate counts in order for a barcode to be considered as lineage representative)
Run 'AggregateCountsHist.R', check the pdf file
Function 'filter' uses the threshold in aggregate counts and outputs filtered versions of all libraries
Function 'trajectories' assembles all barcode filtered counts per assay and outputs a single file (time trajectories with just raw counts): these are not normalized counts, use as input in fitness estimation algorithm
Function 'norm_trajectories' are also the timecourse trajectories but normalized based on the total counts. I plot these to decide on which timepoints to use for fitness estimation
