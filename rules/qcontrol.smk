import os

if not config['dataset']['nanopore']:

    # Assess the quality of demultiplexed Illumina FASTQ files using FastQC.
    # FastQC analyzes sequencing quality metrics, including per-base quality,
    # GC content, sequence duplication, and potential adapter contamination.
    # Each input FASTQ file produces an HTML report and a ZIP archive
    # containing the detailed quality control results.
    rule fastqc:
        input:
            os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_R{read}.fastq.gz")
        output:
            os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.html"),
            os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.zip")
        params:
            config["general"]["output_dir"]
        threads: 1
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/qcontrol/fastqc/{sample}_{unit}_R{read}.log")
        conda:
            "../envs/preprocessing/qcontrol.yaml"
        shell:
            "fastqc -o {params}/qc/ {input} -t {threads} > {log} 2>&1"

    # Combine the individual FastQC results into a single MultiQC report.
    # MultiQC collects the generated FastQC ZIP archives and summarizes
    # their quality metrics to facilitate comparisons across samples.
    # The rule creates a combined HTML report without modifying the reads.
    rule multiqc:
        input:
            expand(os.path.join(config["general"]["output_dir"],"qc/{unit.sample}_{unit.unit}_R{read}_fastqc.zip"),
            unit=units.reset_index().itertuples(), read=reads)
        output:
            os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html")
        params:
            config["general"]["output_dir"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/qcontrol/multiqc.log")
        conda:
            "../envs/preprocessing/qcontrol.yaml"
        shell:
            "multiqc {params}/qc/ -o {params}/qc/ > {log} 2>&1"


else:

    # Assess the quality of input Nanopore FASTQ files using FastQC.
    # FastQC analyzes sequencing quality metrics, including per-base quality,
    # GC content, sequence duplication, and potential adapter contamination.
    # Each input FASTQ file produces an HTML report and a ZIP archive
    # containing the detailed quality control results.
    rule fastqc:
        input:
            os.path.join(config["general"]["filename"],"{sample}_{unit}_R{read}.fastq.gz")
        output:
            os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.html"),
            os.path.join(config["general"]["output_dir"],"qc/{sample}_{unit}_R{read}_fastqc.zip")
        params:
            config["general"]["output_dir"]
        threads: 1
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/qcontrol/fastqc/{sample}_{unit}_R{read}.log")
        conda:
            "../envs/preprocessing/qcontrol.yaml"
        shell:
            "fastqc -o {params}/qc/ {input} -t {threads} > {log} 2>&1"


    # Combine the individual FastQC results into a single MultiQC report.
    # MultiQC collects the generated FastQC ZIP archives and summarizes
    # their quality metrics to facilitate comparisons across samples.
    # The rule creates a combined HTML report without modifying the reads.
    rule multiqc:
        input:
            expand(os.path.join(config["general"]["output_dir"],"qc/{unit.sample}_{unit.unit}_R{read}_fastqc.zip"),
                unit=units.reset_index().itertuples(), read=reads)
        output:
            os.path.join(config["general"]["output_dir"],"qc/multiqc_report.html")
        params:
            config["general"]["output_dir"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/qcontrol/multiqc.log")
        conda:
            "../envs/preprocessing/qcontrol.yaml"
        shell:
            "multiqc {params}/qc/ -o {params}/qc/ > {log} 2>&1"
