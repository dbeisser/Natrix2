if config['dataset']['nanopore']:
    if config['nanopore']['pychopper']:
    ## fastq to fasta
        rule fastq2fasta:
            input:
                expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_merged/{{sample}}_{{unit}}_R{read}.fastq"),read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
            conda:
                "../envs/seqtk.yaml"
            shell:
                "seqtk seq -a {input} > {output}"

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
                " cd-hit-est -i {input} -o {output} -l {params.length} -c 0.8 -d 0 -M {params.memory} -T {threads}"

        ## 1st round of alignment with minimap on reads
        rule minimap_align:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 2nd round of alignment with minimap on reads
        rule minimap_align_2:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 3rd round of alignment with minimap on reads
        rule minimap_align_3:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 4th round of alignment with minimap on reads
        rule minimap_align_4:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.tmp"),read=reads),
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## medaka polishing
        rule medaka_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"), read=reads)
            output:
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                tmp = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)
            conda: "../envs/medaka.yaml"
            params:
                prefix="consensus",
                out_dir=expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}"), read=reads),
                memory=config["general"]["memory"]
            shell:
                """
                sed "s/[>].*[^|]|/>/" {input[0]} > {output.tmp};
                export CUDA_VISIBLE_DEVICES=''
                medaka_consensus -i {output.tmp} -d {input[1]} -o {params.out_dir} -t {threads} -b {params.memory}
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
                "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

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

    else:
        ## fastq to fasta
        rule fastq2fasta:
            input:
                expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads)
            conda:
                "../envs/seqtk.yaml"
            shell:
                "seqtk seq -a {input} > {output}"

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
                " cd-hit-est -i {input} -o {output} -l {params.length} -c 0.8 -d 0 -M {params.memory} -T {threads}"


         ## 1st round of alignment with minimap on reads
        rule minimap_align:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_1.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 2nd round of alignment with minimap on reads
        rule minimap_align_2:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_1.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_2.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 3rd round of alignment with minimap on reads
        rule minimap_align_3:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_2.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_3.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.tmp"),read=reads),
                final = temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads))
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## 4th round of alignment with minimap on reads
        rule minimap_align_4:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_3.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                temp(expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align_4.sam"),read=reads))
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
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.tmp"),read=reads),
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"),read=reads)
            threads: config["general"]["cores"]
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t {threads} {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## medaka polishing
        rule medaka_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon_4.fasta"), read=reads)
            output:
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                tmp = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)
            threads: config["general"]["cores"]
            params:
                prefix="consensus",
                out_dir=expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}"), read=reads),
                memory=config["general"]["memory"]
            conda: "../envs/medaka.yaml"
            shell:
                """
                sed "s/[>].*[^|]|/>/" {input[0]} > {output.tmp};
                export CUDA_VISIBLE_DEVICES=''
                medaka_consensus -i {output.tmp} -d {input[1]} -o {params.out_dir} -t {threads} -b {params.memory}
                """

        ## aligh medaka consesus with raw reads to get number of reads for each consensus
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
                "minimap2 -ax map-ont -t {threads} {input[0]}  {input[1]} > {output}"

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
