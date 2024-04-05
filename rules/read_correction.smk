import os

if config['dataset']['nanopore']:
    ## fastq to fasta
    rule fastq2fasta:
        input:
            pychopper_merged = expand(os.path.join(config["general"]["output_dir"], "pychopper/output/merged/{{sample}}_{{unit}}_R{read}.fastq"), read=reads) if config['nanopore']['pychopper'] else expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"),read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        conda:
            "../envs/seqtk.yaml"
        shell:
            """
            seqtk seq -a {input} > {output.tmp};
            sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
            """
                            
    ## cdhit
    rule cd_hit:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads)
        threads: config["general"]["cores"]
        params:
            memory=config["general"]["memory"],
            length=config["nanopore"]["min_length"]
        conda: "../envs/read_correction.yaml"
        shell:
            "cd-hit-est -i {input} -o {output} -l {params.length} -c 0.8 -d 0 -M {params.memory} -T {threads}"

    ## 1st round of alignment with minimap on reads
    rule minimap_align:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    ## 1st round of polishing with racon
    rule racon_polishing:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"), read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"), read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]/>rp1;/" {output.tmp} > {output.final}
            """

    ## 2nd round of alignment with minimap on reads
    rule minimap_align_2:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    ## 2nd round of polishing with racon
    rule racon_polishing_2:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp1/>rp2rp1/" {output.tmp} > {output.final}
            """

    ## 3rd round of alignment with minimap on reads
    rule minimap_align_3:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    ## 3rd round of polishing with racon
    rule racon_polishing_3:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp2rp1/>rp3rp2rp1/" {output.tmp} > {output.final}
            """

    ## 4th round of alignment with minimap on reads
    rule minimap_align_4:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    ## 4th round of polishing with racon
    rule racon_polishing_4:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp3rp2rp1/>rp4rp3rp2rp1/" {output.tmp} > {output.final}
            """

    ## 5th round of alignment with minimap on reads
    rule minimap_align_5:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_5.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

    ## 5th round of polishing with racon
    rule racon_polishing_5:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_5.sam"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
        output:
            tmp = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_5.tmp"),read=reads)),
            final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_5.fasta"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
                racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                sed "s/[>]rp4rp3rp2rp1/>rp5rp4rp3rp2rp1/" {output.tmp} > {output.final}
            """

    ## medaka polishing
    rule medaka_polishing:
        input:
            fasta_file=expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
            racon_file=expand(os.path.join(config["general"]["output_dir"], "read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_{racon}.fasta"), read=reads, racon=config['nanopore']['racon'])
        params:
            out_dir=expand(os.path.join(config["general"]["output_dir"], "read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp"), read=reads)
        output:
            temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp/consensus.fasta"), read=reads))
        threads: config["general"]["cores"]
        conda: "../envs/medaka.yaml"
        shell:
            """
            export CUDA_VISIBLE_DEVICES=''
            medaka_consensus -i {input.fasta_file} -d {input.racon_file} -o {params.out_dir} -t {threads}
            """

    ## remove double header
    rule rm_racon_header:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/temp/consensus.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads)
        shell:
            """
            sed "s/[>].*[^ ] />/" {input} > {output}
            """

    ## align medaka consensus with raw reads to get number of reads for each consensus
    rule minimap_medaka:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
            expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
        threads: config["general"]["cores"]
        conda:
            "../envs/read_correction.yaml"
        shell:
            """
            minimap2 -ax map-ont -t {threads} {input[0]} {input[1]} > {output}
            """

    rule counts_minimap:
        input:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}/counts.txt"),read=reads),
            expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}/rep_consensus.fasta"),read=reads)
        conda: "../envs/samtools.yaml"
        script:
            "../scripts/counts_consensus_repeat.py"
