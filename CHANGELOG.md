# Change Log
All significant changes to this project are documented here.

## 2024-07-01
- File: chim_rm.smk; The parameter 'threads' has been set to 1 in the rule 'vsearch_chim'.

## 2024-06-27
- Added: Tool fastq_inspector.py
- Added seqkit to natrix.yaml
- README: New section (Sequence Count); Contents added;
- Adjustments in all existing configuration files

## 2024-06-04
- Added parameter: vsearch_id in configfiles: Illumina.yaml, Illumina_swarm.yaml;
- Added dependencies: tool pigz to files: natrix.yaml, blast.yaml;
- Replaced tool: gunzip with pigz; pigz uses multiple CPU cores
- Bug fixed: wget command in file: blast.smk; rule: make_silva_db
- Organized files and folders in the root directory
- Updated tool versions: mothur, fastqc, multiqc
- Added: Error message for ASV error in Snakefile

## 2023-09-22
- Added PR2, UNITE Databases
- Added mumu for post clustering
- Added MOTHUR 
- Added Nanopore workflow

## 2021-03-12
- Cluster split samples together in DADA2
- Updated all packages except SWARM

## 2021-02-19
- Updated Snakemake version

## 2021-01-27
- Changed RAM usage of DADA2 by running clustering on single samples, while still estimating errors on all samples.

## 2020-02-28
- Updated to NCBI version 5 BLAST databases (dbV5) which require different files and BLAST+ > 2.9.0
- Changed SILVA rules accordingly to work with new BLAST
- Older NCBI database versions will no longer be supported (v4 soon depricated at NCBI)

## 2020-03-23
- Added support for amplicon sequence variants using the DADA2 algorithm

## 2020-04-01
- Added more logging information