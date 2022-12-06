# Understanding genetic variation in diverse human populations through transcriptome and genome sequencing

## Retrieving the data locally:

(1) Download the data file ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz from https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/

(2) Download the file GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz from https://zenodo.org/record/6998955#.Y46Ix-zMLzc

(3) Download the file ALL_1000G_phase1integrated_v3.sample from Tiffany (Our capstone domain professor)

Please place these files in data/raw_data

## Running the project

* To install the dependencies launch the dockerhub image using the dockerhub-id specified in submissions.json
* Running python run.py test will output some of the plots included in our report as well as linear regression summaries from our eQTL analyses, and Chromosome, Position, P-Value tables that we use through LocusZoom to create Manhattan plots.
https://my.locuszoom.org/gwas/955311/region/?chrom=22&start=16677262&end=17077262
https://my.locuszoom.org/gwas/466906/region/?chrom=22&start=16705432&end=17105432
https://my.locuszoom.org/gwas/798344/region/?chrom=22&start=16874656&end=17274656
* WARNING: When loading and looping through the VCF you may run out of memory. The VCF contains many dataframe chunks that get filtered to a subset we perform the analysis on. Our test data is a subset of this data used for running Linear Regressions and plotting if you run into issues from the raw data. (Only pertinent for 'raw' run)

## Building the project stages using CMD line arguments 'test' and 'raw'
* Running the command 'python3 run.py raw' will run code and produce outputs from the raw files contained in our data folder under raw_data
* Running the command 'python3 run.py test' will run code and produce outputs using the test data contained in our data folder under raw_data

## Outputs
All outputs from the 'test' and 'raw' runs are located in the output folder.

# DSC180A_Q1Project
