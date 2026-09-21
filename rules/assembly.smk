import os

if not config['dataset']['nanopore']:
    def get_fastq(wildcards):
        if not is_single_end(wildcards.sample, wildcards.unit):
            return expand(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{group}.fastq"),
                            group=[1,2], **wildcards)
        return os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_1.fastq").format(**wildcards)

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
            "seqtk seq -a {input} > {output}"
