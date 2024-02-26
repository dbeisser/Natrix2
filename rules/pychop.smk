import os
if config['dataset']['nanopore'] and config['nanopore']['pychopper']:
    rule define_pychop_primer:
        output:
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        params:
            primer = config["general"]["primertable"]
        shell:
                """
                forward_primer=$(sed -n "2p" {params.primer} | cut -d, -f4);
                reverse_primer=$(sed -n "2p" {params.primer} | cut -d, -f7);
                printf ">SSP\n$forward_primer\n>VNP\n$reverse_primer" > {output.primers};
                echo "+:SSP,-VNP|-:VNP,-SSP" > {output.primers_config}
                """

    rule pychop:
        input:
            fastq=expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        output:
            out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/{{sample}}_{{unit}}_R{read}.tmp"), read=reads),
            pdf_out=expand(os.path.join(config["general"]["output_dir"],"pychopper/reports/{{sample}}_{{unit}}_R{read}.pdf"), read=reads),
            unclass_fastq= expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/{{sample}}_{{unit}}_R{read}.tmp"), read=reads),
            rescue_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/rescued/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            length_out = expand(os.path.join(config["general"]["output_dir"], "pychopper/unclassified/{{sample}}_{{unit}}_R{read}_length_out.tmp"), read=reads)
        threads: config["general"]["cores"]
        params:
            qual = config["nanopore"]["pychopqual"],
            length = config["nanopore"]["min_length"]
        conda:
            "../envs/pychopper.yaml"
        shell:
                """
                pychopper -m edlib -b {input.primers} -Q {params.qual} -z {params.length} -l {output.length_out} -c {input.primers_config} -r {output.pdf_out} -u {output.unclass_fastq} -w {output.rescue_fastq} -t {threads} {input.fastq} {output.out_fastq};
                """

    rule pychopper_rescue:
        input:
            unclass_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/{{sample}}_{{unit}}_R{read}.tmp"), read=reads),
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        output:
            unclass_out_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_unclass/output/{{sample}}_{{unit}}_R{read}.tmp"), read=reads),
            unclass_unclass_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_unclass/unclassified/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads),
            unclass_rescue_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_unclass/rescued/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            unclass_pdf=expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_unclass/reports/{{sample}}_{{unit}}_R{read}.pdf"),  read=reads),
            length_out = expand(os.path.join(config["general"]["output_dir"], "pychopper/unclassified/{{sample}}_{{unit}}_R{read}_length_out.tmp"), read=reads)
        threads: config["general"]["cores"]
        params:
            qual =config["nanopore"]["pychopqual"],
            length = config["nanopore"]["min_length"]
        conda:
            "../envs/pychopper.yaml"
        shell:
            """
            pychopper -m edlib -x rescue -b {input.primers} -Q {params.qual} -z {params.length} -l {output.length_out} -c {input.primers_config} -r {output.unclass_pdf} -u {output.unclass_unclass_fastq} -w {output.unclass_rescue_fastq} -t {threads} {input.unclass_fastq} {output.unclass_out_fastq};
            """

    rule merge_pychopper:
        input:
            out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/{{sample}}_{{unit}}_R{read}.tmp"),  read=reads),
            unclass_out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/pychopper_unclass/output/{{sample}}_{{unit}}_R{read}.tmp"),  read=reads)
        output:
            pychopper_merged = expand(os.path.join(config["general"]["output_dir"], "pychopper/pychopper_merged/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        shell:
            """
            cat {input.out_fastq} {input.unclass_out_fastq} > {output.pychopper_merged}
            """



