import os

# Cluster nucleotide sequences with CD-HIT-EST using the configured identity,
# length-overlap, memory, and thread settings; retain its cluster assignments.
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
        # Run CD-HIT-EST with the configured identity, threads, length coverage and memory limits.
        # Write representative sequences and cluster assignments to the declared outputs.
        # Log stdout and stderr of CD-HIT-EST to the log file and terminal.
        "cd-hit-est -i {input} -o {output[0]} -c {params.id_percent} -T"
        " {threads} -s {params.length_cutoff} -M {params.memory} -sc 1 -d 0"
        " 2>&1 | tee {log:q}"

# Generate the dereplicated FASTA from CD-HIT results and source sequences
# using the configured representative-selection method.
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
