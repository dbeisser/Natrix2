import os

# Process Illumina FASTQ reads and generate assembled sequence files.
# The workflow prepares primer settings, quality-filters reads with PRINSEQ,
# processes reads with Cutadapt, assembles reads using the configured
# assembly settings, and converts the resulting FASTQ files to FASTA
# for subsequent sequence processing.

# [Dataset-specific workflow: Illumina data]
if not config['dataset']['nanopore']:

    # Return the expected demultiplexed FASTQ input path(s) for a sample/unit.
    # Return both numbered reads (1 and 2) for paired-end units; otherwise return only read 1.
    def get_fastq(wildcards):
        if not is_single_end(wildcards.sample, wildcards.unit):
            return expand(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{group}.fastq"),
                            group=[1,2], **wildcards)
        return os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_1.fastq").format(**wildcards)

    # Prepare the primer table used by downstream read-processing rules.
    # The primer information from the user-provided primer table is combined
    # with the configured settings (paired-end mode, primer offset, and whether
    # barcodes or all primers have already been removed).
    rule define_primer:
        input:
            primer_table=config["general"]["primertable"]
        output:
            os.path.join(config["general"]["output_dir"],"primer_table.csv")
        params:
           paired_end=config["merge"]["paired_End"],
           offset=config["qc"]["primer_offset"],
           bar_removed=config["qc"]["barcode_removed"],
           all_removed=config["qc"]["all_primer"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/assembly/define_primer.log")
        conda:
            "../envs/preprocessing/define_primer.yaml"
        script:
            "../scripts/preprocessing/define_primer.py"

    # Process demultiplexed reads with the PRINSEQ helper script.
    # Pass the configured quality threshold to the helper script to
    # produce FASTQ files for subsequent processing.
    # The filtered reads are temporary and deleted once no longer needed.
    rule prinseq:
        input:
            sample=get_fastq
        output:
            temp(expand(
            os.path.join(config["general"]["output_dir"],"assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq"),
            read=reads))
        params:
            config["qc"]["mq"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/assembly/prinseq/{sample}_{unit}.log")
        conda:
            "../envs/preprocessing/prinseq.yaml"
        script:
            "../scripts/preprocessing/prinseq.py"

    # Process quality-filtered reads using the Cutadapt helper script
    # and the generated primer table.
    # Pass the configured minimum and maximum lengths to the script.
    # The trimmed reads are temporary and deleted once no longer needed.
    rule cutadapt:
        input:
            expand(
            os.path.join(config["general"]["output_dir"],"assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq"),
            read=reads),
            primer_t=os.path.join(config["general"]["output_dir"],"primer_table.csv")
        output:
            temp(expand(
            os.path.join(config["general"]["output_dir"],"assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}_cut.fastq"),
            read=reads))
        params:
            paired_end=config["merge"]["paired_End"],
            bar_removed=config["qc"]["barcode_removed"],
            prim_rm=config["qc"]["all_primer"],
            minlen=config["qc"]["minlen"],
            maxlen=config["qc"]["maxlen"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/assembly/cutadapt/{sample}_{unit}.log")
        conda:
            "../envs/preprocessing/cutadapt.yaml"
        script:
            "../scripts/preprocessing/cutadapt.py"

    # Generate assembled FASTQ reads using the assembly helper script.
    # Pass paired-end, overlap, quality, and length settings to the script;
    # the generated primer table is also provided as input.
    # The behavior for single-end units is determined by the helper script.
    # The result is one assembled FASTQ file per sample/unit.
    rule assembly:
        input:
            expand(
            os.path.join(config["general"]["output_dir"],"assembly/{{sample}}_{{unit}}/{{sample}}_{{unit}}_{read}.fastq"),
            read=reads),
            primer_t=os.path.join(config["general"]["output_dir"],"primer_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq")
        params:
            paired_end=config["merge"]["paired_End"],
            threshold=config["qc"]["threshold"],
            minoverlap=config["qc"]["minoverlap"],
            minlen=config["qc"]["minlen"],
            maxlen=config["qc"]["maxlen"],
            minqual=config["qc"]["minqual"],
            prim_rm=config["qc"]["all_primer"]
        threads: 20
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/assembly/assembly/{sample}_{unit}.log")
        conda:
            "../envs/analysis/assembly.yaml"
        script:
            "../scripts/analysis/assembly.py"

    # Convert assembled FASTQ to FASTA using seqtk.
    # The FASTA output contains sequences without FASTQ quality scores.
    # It is provided for downstream sequence processing.
    # The FASTA file is temporary and deleted once no longer needed.
    rule copy_to_fasta:
        input:
            os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_assembled.fastq")
        output:
            temp(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.fasta"))
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/assembly/copy_to_fasta/{sample}_{unit}.log")
        conda:
            "../envs/utilities/seqtk.yaml"
        shell:
            """
                # Convert FASTQ to FASTA (stderr to log, sequences to FASTA file)
                seqtk seq -a {input:q} 2>&1 > {output:q} | tee {log:q}
            """
