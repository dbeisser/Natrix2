import os
if config['dataset']['nanopore']:
    rule pigz_unzip:
        input:
            os.path.join(config["general"]["filename"],"{sample}_{unit}_R{read}.fastq.gz")
        output:
            temp(os.path.join(config["general"]["output_dir"],"fastq/{sample}_{unit}_{read}.tmp"))
        shell:
            "pigz -d -c {input} > {output}"
            
    rule check_fastq_format:
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

    rule quality_filtering:
        input:
            expand(os.path.join(config["general"]["output_dir"],"fastq/{{sample}}_{{unit}}_{read}.tmp"), read=reads)
        output:
            expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        threads: config["general"]["cores"]
        params:
            qual = config["nanopore"]["quality_filt"],
            minlength = config["nanopore"]["min_length"],
            maxlength = config["nanopore"]["max_length"],
            headcrop = config["nanopore"]["head_trim"],
            tailcrop = config["nanopore"]["tail_trim"],
            threads = config["general"]["cores"]
        conda:
            "../envs/quality_filtering.yaml"
        shell:
                """
                cat {input} | chopper --quality {params.qual} --minlength {params.minlength} --maxlength {params.maxlength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} --threads {threads} > {output}
                """



