# DBNascent_Analysis
This repository contains notebooks used to generate figures and intermediate files for the project. These are meta-analyses performed using data in `DBNascent`. 

# Linked Repositories
## Construction of the MySQL database `DBNascent` 

Repository for "building, updating, and querying DBNascent".

https://github.com/Dowell-Lab/DBNascent-build

## Pre-processing
### `Nascent-Flow`

The nextflow pipeline used to preprocess raw files 
* Read quality assessment
* Preprocess reads based on quality
* Read mapping to reference genomes

https://github.com/Dowell-Lab/Nascent-Flow

### `Bidirectional-Flow`

A nextflow pipeline for bidirectional transcript detection 
* Generate input files to bidirectional detection tools
* Detect bidirectional transcripts using:
  * dREG
  * Tfit
* Count reads over gene and bidirectional transcripts

https://github.com/Dowell-Lab/Bidirectional-Flow

## Merge bidirectional calls to create consensus annotations.

Using muMerge and bedtools utilities, bidirectional transcript annotations across all human and mouse samples were merged.

https://github.com/Dowell-Lab/bidirectionals_merged 

## Identifying gene & bidirectional pairs

Correlation based method to link bidirectional transcripts to putative targets.

https://github.com/Dowell-Lab/bidir_gene_pairs
