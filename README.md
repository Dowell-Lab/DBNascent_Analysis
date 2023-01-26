# DBNascent_Analysis
Contains meta analyses performed using data in `DBNascent`

# Linked Repositories

## Processing
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

## Construction of `DBNascent` 

https://github.com/Dowell-Lab/DBNascent-build

## Identifying gene & bidirectional pairs

https://github.com/Dowell-Lab/bidir_gene_pairs
