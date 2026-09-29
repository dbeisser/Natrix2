import os

rule cdhit:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.fasta") if not config['dataset']['nanopore'] else expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}/rep_consensus.fasta"),read=reads)
    output:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta"),
        temp(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr"))
    params:
        id_percent=config["derep"]["clustering"],
        length_cutoff=config["derep"]["length_overlap"],
        memory=config["general"]["memory"]
    threads: config["general"]["cores"]
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/dereplication/cdhit/{sample}_{unit}.log")
    conda:
        "../envs/analysis/dereplication.yaml"
    shell:
        "cd-hit-est -i {input} -o {output[0]} -c {params.id_percent} -T"
        " {threads} -s {params.length_cutoff} -M {params.memory} -sc 1 -d 0"

rule cluster_sorting:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta"),
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_cdhit.fasta.clstr"),
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.fasta") if not config['dataset']['nanopore'] else expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"),read=reads)
    output:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta")
    params:
        repr=config["derep"]["representative"],
        length_cutoff=config["derep"]["length_overlap"]
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/dereplication/cluster_sorting/{sample}_{unit}.log")
    conda:
        "../envs/analysis/dereplication.yaml"
    script:
        "../scripts/analysis/dereplication.py"
