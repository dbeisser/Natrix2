if config['dataset']['nanopore']:
    if config['dataset']['pychopper']:
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
            conda: "../envs/read_correction.yaml"
            shell:
                " cd-hit-est -i {input} -o {output} -c 0.9 -d 0 -M 0 -T 0"

        ## alighnment with minimap on reads
        rule minimap_align:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
            conda:
                 "../envs/read_correction.yaml"
            shell:
                "minimap2 -ax map-ont -t 20 {input[0]}  {input[1]} > {output}"

        ## polishing with racon
        rule racon_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"), read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"), read=reads)
            output:
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.tmp"),read=reads),
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.fasta"),read=reads)
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t 20 {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## medaka polishing
        rule medaka_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.fasta"), read=reads)
            params:
                prefix="consensus",
                out_dir=expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}"), read=reads)
            output:
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                tmp = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)
            conda: "../envs/medaka.yaml"
            shell:
                """
                sed "s/[>].*[^|]|/>/" {input[0]} > {output.tmp};
                export CUDA_VISIBLE_DEVICES=''
                medaka_consensus -i {output.tmp} -d {input[1]} -o {params.out_dir} -t 20
                """
    ## aligh medaka consesus with raw reads to get number of reads for each consensus

        rule minimap_medaka:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
            conda:
                 "../envs/read_correction.yaml"
            shell:
                "minimap2 -ax map-ont -t 20 {input[0]}  {input[1]} > {output}"

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
        ## gunzip fastq
        rule gunzip:
            input:
                os.path.join(config["general"]["filename"],"{sample}_{unit}_{read}.fastq.gz")
            output:
                temp(os.path.join(config["general"]["output_dir"],"fastq/{sample}_{unit}_{read}.tmp"))
            shell:
                "gunzip -c {input} > {output}"
        
        ## check format        
        rule check_file_format:
            input:
                os.path.join(config["general"]["output_dir"],"fastq/{sample}_{unit}_{read}.tmp")
            output:
                temp(os.path.join(config["general"]["output_dir"],"fastq/{sample}_{unit}_{read}.fastq"))
            shell:
                """
                    if find {input} -not -type d -exec file '{{}}' ';' | grep CRLF
                    then
                        sed 's/\r$//' {input} > {output}
                    else
                        mv {input} {output}
                    fi
                """

        ## fastq to fasta
        rule fastq2fasta:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fastq/{{sample}}_{{unit}}_R{read}.fastq"),read=reads)
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
            conda: "../envs/read_correction.yaml"
            shell:
                " cd-hit-est -i {input} -o {output} -c 0.9 -d 0 -M 0 -T 0"

        ## alighnment with minimap on reads
        rule minimap_align:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
            conda:
                 "../envs/read_correction.yaml"
            shell:
                "minimap2 -ax map-ont -t 20 {input[0]}  {input[1]} > {output}"

        ## polishing with racon
        rule racon_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/minimap/{{sample}}_{{unit}}_R{read}_align.sam"), read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/cd_hit/{{sample}}_{{unit}}_R{read}_rep.fasta"), read=reads)
            output:
                tmp = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.tmp"),read=reads),
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.fasta"),read=reads)
            conda:
                "../envs/read_correction.yaml"
            shell:
                """
                    racon  -m 8 -x -6 -g -8 -w 500 -t 20 {input[0]} {input[1]} {input[2]} > {output.tmp};
                    sed "s/[>].*[^|]|/>/" {output.tmp} > {output.final}
                """

        ## medaka polishing
        rule medaka_polishing:
            input:
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"),read=reads),
                expand(os.path.join(config["general"]["output_dir"],"read_correction/racon/{{sample}}_{{unit}}_R{read}_racon.fasta"), read=reads)
            params:
                prefix="consensus",
                out_dir=expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}"), read=reads)
            output:
                final = expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                tmp = expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.tmp"),read=reads)
            conda: "../envs/medaka.yaml"
            shell:
                """
                sed "s/[>].*[^|]|/>/" {input[0]} > {output.tmp};
                export CUDA_VISIBLE_DEVICES=''
                medaka_consensus -i {output.tmp} -d {input[1]} -o {params.out_dir} -t 20
                """
         ## aligh medaka consesus with raw reads to get number of reads for each consensus

        rule minimap_medaka:
            input:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/medaka/{{sample}}_{{unit}}_R{read}/consensus.fasta"), read=reads),
                expand(os.path.join(config["general"]["output_dir"],"fasta/{{sample}}_{{unit}}_R{read}.fasta"), read=reads)
            output:
                expand(os.path.join(config["general"]["output_dir"],"read_correction/counts_mapping/{{sample}}_{{unit}}_R{read}_align.sam"),read=reads)
            conda:
                "../envs/read_correction.yaml"
            shell:
                "minimap2 -ax map-ont -t 20 {input[0]}  {input[1]} > {output}"

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
