#### Natrix2; Bioinformatics
# Changelog
This document records all notable changes to Natrix2 project.

## [2025-05-30]
- Update CHANGELOG.md and README.md
- Update Dockerfile to fix build issue
- Update version: seqkit=2.10.0
- Release Natrix2 v2.0.0: major update
- Rename tool fastq_inspector.py to nseqc.py  
- Improve nseqc.py with code documentation

## [2025-05-13]
- Fix OTU merge bug in merge_results.py

## [2025-05-04]
- Remove demultiplexing option for illumina datatype
- Update config files
- Add syntax comments to Snakefile
- Rename folder tools/ to natrixlib/
- Update README.md

## [2025-04-22]
- Update pipeline.sh; improve overall script structure
- Add screen dependency to environment in natrix.yaml
- Fix OTU merging issue in Mothur with SILVA; add rule filter_unclassified
- Add script filter_unclassified.py for taxonomy cleanup
- Update README.md and CHANGELOG.md
- Create folder (documentation/manuals) for user manuals

## [2025-04-08]
- Update main config structure
- Create folder (config_presets/) for pre-set configs
- Move and update test configs to config_presets/
- Remove outdated test configs
- Update Docker dummy files
- Add Mothur CutOff parameter to config

## [2024-12-27]
- Change sample names in input_data/Nanopore_data folder

## [2024-12-19]
- Enforce valid sample names with clear error message
- Revise DADA2 rule for efficient resource utilization
- Protect output files in the DADA2 rule if an error occurs
- Remove redundant loop in DADA2 script

## [2024-12-05]
- Fix clustering logic in classify.smk for ASV
- Update version: bioconductor-dada2=1.30.0

## [2024-11-20]
- Improve velocity with pigz in applicable rules
- Update database versions: SILVA=138.2, PR2=5.0.0
- Adjust MultiQC version to prevent issues
- Adapt configuration files for new database version

## [2024-11-04]
- Add docker_manual.pdf and docker_manual.tex
- Update docker-compose.yaml
- Update README.md

## [2024-10-24]
- Add test_docker.yaml for Docker container testing
- Update docker_pipeline.sh script

## [2024-10-16]
- Update version: multiqc=1.25.1
- Update README.md
- Update docker_pipeline.sh, docker-compose.yaml, natrix.yaml
- Build new Docker Image; natrix2:latest
- Update Natrix2 version: v1.1.2

## [2024-10-04]
- Add file verification: check for uncompressed files
- Add folder to .gitignore
- README: Add new section: Use dev-branch
- Update README.md
- Update natrix.yaml

## [2024-09-30]
- Update filename for ./pipeline.sh to fix execution

## [2024-09-24]
- Add agreement note to the download link of the UNITE Database
- Add new Docker Image release: Natrix2: v1.1.1
- Update README.md
- Update Docker Image, Dockerfile and docker-compose
- Update SILVA-Database Download

## [2024-08-19]
- Update README.md

## [2024-08-12]
- Update UNITE Database; Release date: 2024-04-04; Version: 10.0

## [2024-07-25]
- Add file: .gitignore

## [2024-07-22]
- Bug fix: Adjust the rule make_silva_db to dynamically download the latest SILVA database
- Add tool: fastq_inspector.py; see section Sequence Count
- Add package: seqkit=2.8.2 to natrix.yaml
- Add package: python=3.8 to filtering.yaml
- Add package: pigz=2.8 to blast.yaml and natrix.yaml
- Replace tool: gunzip with pigz in quality_filt.smk, demultiplexing.smk, blast.smk
- Add folders: input_data, primer_table, dag_plots, tools
- README: Add new sections: Sequence Count and Table of Contents
- Organize files and folders in the root directory
- Adjust all existing configuration files
- Update tool versions: fastqc=0.12.1, multiqc=1.23

## [2024-07-18]
- Add checks to ensure seq_rep: ASV is not used with incompatible options in create_dataframe.py
- Fix issue in error message handling in create_dataframe.py

## [2024-07-02]
- Add channels for package availability

## [2024-07-01]
- Fix issues in read_assembly.smk and assembly.py
- Set threads to 1 in vsearch_chim rule in chim_rm.smk

## [2024-06-24]
- Update cutadapt.py: Change the order of barcode and primer pattern processing

## [2023-09-22]
- Add PR2, UNITE Databases
- Add mumu for post clustering
- Add MOTHUR
- Add Nanopore workflow

## [2021-03-12]
- Cluster split samples together in DADA2
- Update all packages except SWARM

## [2021-02-19]
- Update Snakemake version

## [2021-01-27]
- Change RAM usage of DADA2 by running clustering on single samples, while still estimating errors on all samples

## [2020-04-01]
- Add more logging information

## [2020-03-23]
- Add support for amplicon sequence variants using the DADA2 algorithm

## [2020-02-28]
- Update to NCBI version 5 BLAST databases (dbV5) which require different files and BLAST+ > 2.9.0
- Change SILVA rules accordingly to work with new BLAST
- Older NCBI database versions will no longer be supported (v4 soon deprecated at NCBI)