# Change Log
All significant changes to this project are documented here.

## [2024-07-18]
- Add checks to ensure seq_rep: ASV is not used with incompatible options in 'create_dataframe.py'
- Fix issue in error message handling in 'create_dataframe.py'

## [2023-09-22]
- Add PR2, UNITE Databases.
- Add mumu for post clustering.
- Add MOTHUR.
- Add Nanopore workflow.

## [2021-03-12]
- Cluster split samples together in DADA2.
- Update all packages except SWARM.

## [2021-02-19]
- Update Snakemake version.

## [2021-01-27]
- Change RAM usage of DADA2 by running clustering on single samples, while still estimating errors on all samples.

## [2020-02-28]
- Update to NCBI version 5 BLAST databases (dbV5) which require different files and BLAST+ > 2.9.0.
- Change SILVA rules accordingly to work with new BLAST.
- Older NCBI database versions will no longer be supported (v4 soon deprecated at NCBI).

## [2020-03-23]
- Add support for amplicon sequence variants using the DADA2 algorithm.

## [2020-04-01]
- Add more logging information.