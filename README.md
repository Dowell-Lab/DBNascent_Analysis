# DBNascent_Analysis
This repository contains notebooks used to generate figures and intermediate files for the project. These are meta-analyses performed using data in `DBNascent`. 

## Structure of repository

- **analysis**: Jupyter notebooks used for analysis and figure generation

- **data**: Smaller data used for analysis (See Zenodo for larger datasets 10.5281/zenodo.10223322)

- **scripts**: Scripts used for analysis

  - `GENIE3_bidir_gene_pairs.R`: SPECs scores python script

  - `GENIE3_bidir_gene_pairs.R`: GENIE3 R script

  - `random_bidirectionals_sampled.R`: Shuffling bidirectionals for relative false positive rates

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
