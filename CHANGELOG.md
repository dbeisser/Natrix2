#### Natrix2; bioinformatics

# Changelog

This document records all notable changes to the **Natrix2** project.  
Use the **dev** branch for the latest updates and features.

---

### [2026-01-27]
- Apply `filter_unclassified` to Nanopore Mothur outputs
- Update `classify.smk`
- Complete test run for nanopore data

---

### [2026-01-01]
- Update README structure and content
- Update CHANGELOG
- Bump MultiQC to v1.33 and fix Conda environment issue
- Refine `filter_unclassified.py` unclassified handling
- Apply `filter_unclassified` to Illumina Mothur outputs
- Fix logical operator in `classify.smk` rule for vsearch
- Reorder channels in `natrix.yaml`
- Fix Mothur `_unclassified` taxonomy repetition
- Release Natrix2 v2.0.1

---

### [2025-12-03]
- Rewrite `merge_mumu_output_2.py` with new merge-key logic
- Add handling of empty/missing mumu OTUs
- Preserve all BLAST OTUs (right-join)
- Generate new OTU identifiers from merge key and size
- Add `unmerged_seqids.csv` output
- Standardize outputs to use `seqid` as index label
- Clean up script and update Snakemake rule

---

### [2025-11-17]
- Fix download issue of the NCBI database for BLAST
- Add packages to `blast.yaml`
- Update BLAST version in `blast.yaml`
- Remove the `screen` package from `natrix.yaml`
- Update and revise `CHANGELOG.md`

---

### [2025-10-12]
- Update `README.md` and revise the Troubleshooting section

---

### [2025-09-26]
- Add new section "Common issues" to `README`
  - Add "Causes of pipeline aborts"
  - Add "Problems with installation"
- Update `.gitignore`

---

### [2025-07-14]
- Update text color in `filename.png` to improve readability
- Add citation section to `README.md`
- Add `CITATION.cff` file for GitHub citation support

---

### [2025-06-13]
- Update `README.md` with the latest usage instructions and options
- Update `docker_manual.tex` and `docker_manual.pdf`

---

### [2025-06-12]
- Add output file `multi_hits_blast_mumu.csv` to log all duplicate OTUs
- Allow `no_filter` as a valid option for `drop_tax_classes` in the config file
- Update merge script `merge_mumu_output_2.py`:
  - Retain only duplicate entries with pident ≥ 90
  - Select the best hit per OTU based on highest pident
  - Export all duplicates to `multi_hits_blast_mumu.csv` for traceability
- Update all pre-set config files
- Update `README.md` with the latest usage instructions and options
- Improve image backgrounds for better visual clarity

---

### [2025-05-30]
- Update `CHANGELOG.md` and `README.md`
- Update `Dockerfile` to fix build issue
- Update version: `seqkit=2.10.0`
- Release Natrix2 v2.0.0 (major update)
- Rename tool `fastq_inspector.py` to `nseqc.py`
- Improve `nseqc.py` with code documentation

---

### [2025-05-13]
- Fix OTU merge bug in `merge_results.py`

---

### [2025-05-04]
- Remove demultiplexing option for Illumina datatype
- Update config files
- Add syntax comments to `Snakefile`
- Rename folder `tools/` to `natrixlib/`
- Update `README.md`

---

### [2025-04-22]
- Update `pipeline.sh` and improve overall script structure
- Add `screen` dependency to `natrix.yaml`
- Fix OTU merging issue in Mothur with SILVA; add rule `filter_unclassified`
- Add script `filter_unclassified.py` for taxonomy cleanup
- Update `README.md` and `CHANGELOG.md`
- Create folder `documentation/manuals` for user manuals

---

### [2025-04-08]
- Update main config structure
- Create folder `config_presets/` for pre-set configs
- Move and update test configs into `config_presets/`
- Remove outdated test configs
- Update Docker dummy files
- Add Mothur CutOff parameter to config

---

### [2024-12-27]
- Change sample names in `input_data/Nanopore_data` folder

---

### [2024-12-19]
- Enforce valid sample names with clear error message
- Revise DADA2 rule for efficient resource utilization
- Protect output files in the DADA2 rule if an error occurs
- Remove redundant loop in DADA2 script

---

### [2024-12-05]
- Fix clustering logic in `classify.smk` for ASV
- Update version: `bioconductor-dada2=1.30.0`

---

### [2024-11-20]
- Improve velocity with `pigz` in applicable rules
- Update database versions: `SILVA=138.2`, `PR2=5.0.0`
- Adjust MultiQC version to prevent issues
- Adapt configuration files for new database versions

---

### [2024-11-04]
- Add `docker_manual.pdf` and `docker_manual.tex`
- Update `docker-compose.yaml`
- Update `README.md`

---

### [2024-10-24]
- Add `test_docker.yaml` for Docker container testing
- Update `docker_pipeline.sh` script

---

### [2024-10-16]
- Update version: `multiqc=1.25.1`
- Update `README.md`
- Update `docker_pipeline.sh`, `docker-compose.yaml`, `natrix.yaml`
- Build new Docker image `natrix2:latest`
- Update Natrix2 version: v1.1.2

---

### [2024-10-04]
- Add file verification to check for uncompressed files
- Add folder to `.gitignore`
- Add new section “Use dev-branch” to `README`
- Update `README.md`
- Update `natrix.yaml`

---

### [2024-09-30]
- Update filename in `pipeline.sh` to fix execution

---

### [2024-09-24]
- Add agreement note to the download link for the UNITE database
- Add new Docker image release Natrix2 v1.1.1
- Update `README.md`
- Update Docker image, `Dockerfile`, and `docker-compose`
- Update SILVA database download

---

### [2024-08-19]
- Update `README.md`

---

### [2024-08-12]
- Update UNITE database (release date: 2024-04-04; version: 10.0)

---

### [2024-07-25]
- Add `.gitignore`

---

### [2024-07-22]
- Fix rule `make_silva_db` to dynamically download the latest SILVA database
- Add tool `fastq_inspector.py` (see section Sequence Count)
- Add package: `seqkit=2.8.2` to `natrix.yaml`
- Add package: `python=3.8` to `filtering.yaml`
- Add package: `pigz=2.8` to `blast.yaml` and `natrix.yaml`
- Replace tool `gunzip` with `pigz` in `quality_filt.smk`, `demultiplexing.smk`, `blast.smk`
- Add folders: `input_data`, `primer_table`, `dag_plots`, `tools`
- Add new sections “Sequence Count” and “Table of Contents” to `README`
- Organize root directory structure
- Adjust all existing configuration files
- Update tool versions: `fastqc=0.12.1`, `multiqc=1.23`

---

### [2024-07-18]
- Add checks to ensure `seq_rep: ASV` is not used with incompatible options in `create_dataframe.py`
- Fix issue in error message handling in `create_dataframe.py`

---

### [2024-07-02]
- Add channels for package availability

---

### [2024-07-01]
- Fix issues in `read_assembly.smk` and `assembly.py`
- Set `threads=1` in `vsearch_chim` rule in `chim_rm.smk`

---

### [2024-06-24]
- Update `cutadapt.py`: Change processing order of barcode and primer patterns

---

### [2023-09-22]
- Add PR2 and UNITE databases
- Add mumu for post-clustering
- Add Mothur
- Add Nanopore workflow

---

### [2021-03-12]
- Cluster split samples together in DADA2
- Update all packages except SWARM

---

### [2021-02-19]
- Update Snakemake version

---

### [2021-01-27]
- Reduce RAM usage in DADA2 by running clustering on single samples while still estimating errors on all samples

---

### [2020-04-01]
- Add additional logging information

---

### [2020-03-23]
- Add support for amplicon sequence variants using the DADA2 algorithm

---

### [2020-02-28]
- Update to NCBI dbV5 BLAST databases requiring BLAST+ > 2.9.0
- Update SILVA rules to work with the new BLAST
- Drop support for older NCBI database versions (v4 soon deprecated)
