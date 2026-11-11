import pandas as pd
import os

# Load sample metadata
# Read sample metadata using (sample, unit) as the two-level index.
# Preserve metadata fields as strings.
units = pd.read_table(os.path.join(config["general"]["output_dir"],config["general"]["units"]), 
    index_col=["sample", "unit"],
    dtype=str)

# Convert index levels to strings
# Normalize index levels to strings for consistent lookup.
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

# Remove read suffix (e.g. R1/R2)
# Remove the final character from the configured read-name suffix.
name_ext = config["merge"]["name_ext"][:-1]

# Check if reads are single-end
# A missing fq2 value identifies a single-end sample/unit.
def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample,unit), "fq2"])

# Define read layout
# Use two read identifiers for paired-end runs; otherwise use the existing scalar.
if config["merge"]["paired_End"]:
    reads = [1,2]
else:
    reads = 1

# Final workflow targets
# The requested outputs determine which rules Snakemake needs to execute.
# Optional outputs become empty lists when their feature is disabled.
rule all:
    input:
        # Require a MultiQC report only when MultiQC is enabled.
        os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html") if config["general"]["multiqc"] else [],
        # Always produce the unfiltered abundance table.
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_table.csv"),
        # Always produce the filtered abundance table.
        os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
        # Save AmpliconDuo results when enabled for split-sample filtering.
        os.path.join(config["general"]["output_dir"],"filtering/figures/AmpliconDuo.RData") if config["merge"]["ampliconduo"] and config["merge"]["filter_method"] == "split_sample" else [],
        # Require the Swarm OTU table for Illumina OTU workflows using Swarm.
        os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and config['clustering']=="swarm" and not config['dataset']['nanopore'] else [],
        # Require VSEARCH OTU abundances in VSEARCH-based OTU mode.
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus_tab.txt") if config["general"]["seq_rep"] == "OTU" and config['clustering']=="vsearch"  else [],
        # Require the VSEARCH table whenever VSEARCH clustering is selected.
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv") if config['clustering']=="vsearch" else [],
        # Create classified tables for all selected MOTHUR databases.
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table.csv"), database=config['classify']['database']) if config['classify']['mothur'] else [],
        # Add MUMU-refined MOTHUR tables when the OTU/MUMU options are enabled.
        expand(os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"), database=config['classify']['database']) if config["general"]["seq_rep"] == "OTU" and  config['classify']['mothur'] and config['postcluster']['mumu'] else [],
        # Create classified tables for the selected BLAST databases.
        expand(os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table.csv"), database=config['blast']['database'].lower()) if config["blast"]["blast"]  else [],
        # Add MUMU-refined BLAST tables when OTU, BLAST and MUMU are enabled.
        expand(os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table_mumu.csv"), database=config['blast']['database'].lower()) if config["general"]["seq_rep"] == "OTU" and  config["blast"]["blast"] and config['postcluster']['mumu'] else [],

# Prefer assembly when multiple rules can generate the same output
# For non-Nanopore runs, prefer assembly if rules overlap in output.
if not config['dataset']['nanopore']:
    ruleorder: assembly > prinseq

# Import workflow rules
# Download and prepare reference databases
include: "rules/databases.smk"
# Demultiplex raw sequencing reads by sample
include: "rules/demultiplexing.smk"
# Generate FastQC and MultiQC quality reports
include: "rules/qcontrol.smk"
# Trim primers and assemble Illumina reads
include: "rules/assembly.smk"
# Detect, orient, and trim full-length Nanopore reads
include: "rules/pychopper.smk"
# Correct Nanopore sequencing errors and generate consensus reads
include: "rules/rcorrection.smk"
# Filter reads by quality and length
include: "rules/qfiltering.smk"
# Dereplicate sequences using CD-HIT
include: "rules/dereplication.smk"
# Detect and remove chimeric sequences
include: "rules/chimera.smk"
# Merge samples and apply abundance filtering
include: "rules/merging.smk"
# Generate OTUs or ASVs using Swarm, DADA2, or VSEARCH
include: "rules/clustering.smk"
# Cluster sequences with VSEARCH
include: "rules/vsearch.smk"
# Merge highly similar OTUs using MUMU
include: "rules/mumu.smk"
# Assign taxonomy using MOTHUR reference databases
include: "rules/classify.smk"
# Assign taxonomy using BLAST reference databases
include: "rules/blast.smk"
