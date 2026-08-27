import os

# Nanopore Read Correction and Polishing Workflow
#
# Pipeline:
#   FASTQ => FASTA => CD-HIT clustering => 5x (minimap2 + racon) => Medaka CNN
#   => Consensus abundance mapping
#
# Design:
#   - CD-HIT collapses redundant reads to provide representative seeds
#   - Each racon round uses the previous consensus as reference, progressively
#     reducing the error rate; minimap2 provides the SAM input racon requires
#   - Medaka (neural network) replaces a 6th racon round for better accuracy
#     on systematic ONT errors (e.g. homopolymers)


if config['dataset']['nanopore']:

    # Convert pychopper/quality-filtered FASTQ => FASTA.
    # seqtk handles format conversion; sed cleans up headers by stripping
    # everything after the last "|" for unambiguous downstream identifiers.
    rule fastq2fasta:
        input:
            # Full-length Nanopore FASTQ reads
            pychopper_merged = expand(os.path.join(config["general"]["output_dir"], "pychopper/output/merged/{{sample}}_{{unit}}_R{read}.fastq"), read=reads) if config['nanopore']['pychopper'] else expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"),read=reads)
        output:
            # Temporary FASTA conversion
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)),
            # FASTA reads for clustering/polishing
            final = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        conda:
            "../envs/utilities/seqtk.yaml"
        shell:
            """
            seqtk seq -a {input} > {output.tmp};
            sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
            """
                            
    # Cluster reads at 80 % identity (CD-HIT-EST) to produce a non-redundant
    # representative set used as the round-1 polishing reference.
    # Params: -l min length, -c identity, -d 0 unlimited header, -M memory, -T threads
    rule cd_hit:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        output:
            # Clustered representative sequences
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads)
        threads: config["general"]["cores"]
        params:
            memory=config["general"]["memory"],
            length=config["nanopore"]["min_length"]
        conda: "../envs/analysis/rcorrection.yaml"
        shell:
            "cd-hit-est -i {input} -o {output} -l {params.length} -c 0.8 -d 0 -M {params.memory} -T {threads}"


# Iterative polishing: rounds 1–5 (minimap2 + racon)
#
# Each round:
# 1. minimap2 maps all raw reads against the current consensus → SAM
# 2. racon uses raw reads + SAM + consensus to compute an improved FASTA
#
# Round 1 reference: CD-HIT representatives
# Rounds 2–5: output of the previous racon round
#
# racon params (fixed across all rounds):
# -m 8 (match); -x -6 (mismatch); -g -8 (gap); -w 500 (window)
#
# sed tracks polishing history in headers (">rp2rp1" = rounds 1+2 done)
# to avoid duplicate identifiers and enable traceability.

    # --- Round 1 ---
    rule minimap_align:
        input:
            # Clustered representative sequences
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for Racon
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    rule racon_polishing:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Minimap2 read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"), read=reads),
            # Clustered representative sequences
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"), read=reads)
        output:
            # Temporary Racon-polished consensus
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.tmp"),read=reads)),
            # Racon-polished consensus FASTA
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]/>rp1;/" {output.tmp} > {output.final}
            """

    # --- Round 2 (ref: _racon_1 => out: _racon_2, header: >rp1 => >rp2rp1) ---
    rule minimap_align_2:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for Racon
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    rule racon_polishing_2:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Minimap2 read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads),
            # Previous Racon-polished consensus
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads)
        output:
            # Temporary Racon-polished consensus
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.tmp"),read=reads)),
            # Racon-polished consensus FASTA
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp1/>rp2rp1/" {output.tmp} > {output.final}
            """

    # --- Round 3 (ref: _racon_2 => out: _racon_3, header: >rp2rp1 => >rp3rp2rp1) ---
    rule minimap_align_3:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for Racon
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    rule racon_polishing_3:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Minimap2 read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads),
            # Previous Racon-polished consensus
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads)
        output:
            # Temporary Racon-polished consensus
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.tmp"),read=reads)),
            # Racon-polished consensus FASTA
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp2rp1/>rp3rp2rp1/" {output.tmp} > {output.final}
            """

    # --- Round 4 (ref: _racon_3 => out: _racon_4, header: >rp3rp2rp1 => >rp4rp3rp2rp1) ---
    rule minimap_align_4:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for Racon
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    rule racon_polishing_4:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Minimap2 read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads),
            # Previous Racon-polished consensus
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads)
        output:
            # Temporary Racon-polished consensus
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.tmp"),read=reads)),
            # Racon-polished consensus FASTA
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp3rp2rp1/>rp4rp3rp2rp1/" {output.tmp} > {output.final}
            """

    # --- Round 5 (ref: _racon_4 => out: _racon_5, header: >rp4rp3rp2rp1 => >rp5rp4rp3rp2rp1) ---
    rule minimap_align_5:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for Racon
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_5.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    rule racon_polishing_5:
        input:
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Minimap2 read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_5.sam"),read=reads),
            # Previous Racon-polished consensus
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
        output:
            # Temporary Racon-polished consensus
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_5.tmp"),read=reads)),
            # Racon-polished consensus FASTA
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_5.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp4rp3rp2rp1/>rp5rp4rp3rp2rp1/" {output.tmp} > {output.final}
            """


    # Neural-network polishing with Medaka as final correction step.
    # Input: all raw reads + racon consensus after N rounds (config: nanopore.racon).
    # CUDA_VISIBLE_DEVICES='' enforces CPU mode for portability.
    # Output: temp consensus.fasta; headers cleaned in rm_racon_header.
    rule medaka_polishing:
        input:
            # Raw FASTA reads
            fasta_file=expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            # Final Racon-polished consensus
            racon_file=expand(os.path.join(config["general"]["output_dir"], "read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_{racon}.fasta"), read=reads, racon=config['nanopore']['racon'])
        params:
            out_dir=expand(os.path.join(config["general"]["output_dir"], "read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp"), read=reads)
        output:
            # Temporary Medaka consensus FASTA
            temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp/consensus.fasta"), read=reads))
        threads: config["general"]["cores"]
        conda: "../envs/analysis/medaka.yaml"
        shell:
            """
            export CUDA_VISIBLE_DEVICES=''
            medaka_consensus -i {input.fasta_file} -d {input.racon_file} -o {params.out_dir} -t {threads}
            """

    # Strip compound racon-history headers from Medaka output (e.g. ">rp5rp4...")
    # to produce clean identifiers for downstream analysis.
    rule rm_racon_header:
        input:
            # Medaka consensus with Racon history headers
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp/consensus.fasta"), read=reads)
        output:
            # Cleaned Medaka consensus FASTA
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads)
        shell:
            """
            sed "s/[>].*[^ ] />/" {input} > {output}
            """

    # Map all raw reads back to Medaka consensus to quantify read support
    # per consensus sequence (used as abundance proxy in counts_minimap).
    rule minimap_medaka:
        input:
            # Final Medaka consensus reference
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
            # Raw FASTA reads
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            # Read-to-consensus alignment for abundance counting
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/analysis/rcorrection.yaml"
        shell:
            """
            minimap2 -ax map-ont -t {threads} {input[0]} {input[1]} > {output}
            """

    # Count raw reads per consensus (abundance proxy) and output a filtered
    # rep_consensus.fasta for downstream taxonomic classification.
    # See scripts/counts_consensus_repeat.py for filtering thresholds.
    rule counts_minimap:
        input:
            # Read-to-consensus alignment
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads),
            # Final Medaka consensus FASTA
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads)
        output:
            # Read counts per consensus sequence
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}/counts.txt"),read=reads),
            # Filtered representative consensus FASTA
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}/rep_consensus.fasta"),read=reads)
        conda: "../envs/utilities/samtools.yaml"
        script:
            "../scripts/analysis/counts_consensus_repeat.py"
