# Change Log
## 2023-09-22
- added PR2, UNITE Databases
- added mumu for post clustering
- added MOTHUR 
- added Nanopore workflow

## 2021-03-12

### Changed
- Cluster split samples together in DADA2
- updated all packages except SWARM

## 2021-02-19

### Changed
- Updated Snakemake version

## 2021-01-27

### Changed
- Changed RAM usage of DADA2 by running clustering on single samples, while still estimating errors on all samples.

## 2020-02-28

### Changed
- updated to NCBI version 5 BLAST databases (dbV5) which require different files and BLAST+ > 2.9.0
- changed SILVA rules accordingly to work with new BLAST
- older NCBI database versions will no longer be supported (v4 soon depricated at NCBI)

## 2020-03-23

### Changed
- Added support for amplicon sequence variants using the DADA2 algorithm

## 2020-04-01

### Changed
- Added more logging information
