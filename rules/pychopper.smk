import os

# Pychopper: Full-length cDNA read detection and rescue
#
# Pipeline:
# Primer definition => Standard classification: primer detection, trimming,
# and orientation (pychop) => Rescue of unclassified reads with relaxed
# alignment (pychopper_rescue) => Merge classified outputs
#
# Pychopper identifies reads containing both expected primers (SSP + VNP),
# trims primer sequences, and re-orients reverse-complement reads into a
# consistent 5' to 3' direction. Reads that fail standard classification are
# re-processed with a rescue algorithm to recover reads with incomplete primer matches.


if config['dataset']['nanopore'] and config['nanopore']['pychopper']:

    # Extract forward (SSP) and reverse (VNP) primers from the primer table
    # and write them to a FASTA file.
    #
    # The config file encodes expected read orientations:
    # "+:SSP,-VNP" = forward strand; "-:VNP,-SSP" = reverse strand
    rule define_pychop_primer:
        output:
            # SSP/VNP primer sequences in FASTA format
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            # Expected primer orientation configuration
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        params:
            # CSV table containing SSP and VNP primer sequences
            primer = config["general"]["primertable"]
        shell:
                """
                forward_primer=$(sed -n "2p" {params.primer} | cut -d, -f4);
                reverse_primer=$(sed -n "2p" {params.primer} | cut -d, -f7);
                printf ">SSP\n$forward_primer\n>VNP\n$reverse_primer" > {output.primers};
                echo "+:SSP,-VNP|-:VNP,-SSP" > {output.primers_config}
                """

    # Standard pychopper run: classifies reads by primer presence, trims
    # primers, orients reads, and filters by quality and minimum length.
    # -m edlib: use edlib aligner for primer detection
    # -Q: minimum quality score (config: pychopqual)
    # -z: minimum read length after trimming (config: min_length)
    # -u: unclassified reads (missing one or both primers) => rescue input
    # -w: lower-confidence reads classified by rescue logic
    # -r: PDF report with classification statistics
    # -l: reads filtered out by length (temp)
    rule pychop:
        input:
            # Quality-filtered Nanopore FASTQ reads
            fastq=expand(os.path.join(config["general"]["output_dir"],"quality_filtering/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            # SSP/VNP primer FASTA for primer detection
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            # Expected primer orientations for read classification
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        output:
            # Classified full-length cDNA reads
            out_fastq = temp(expand(os.path.join(config["general"]["output_dir"],"pychopper/output/normal/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)),
            # Pychopper classification report
            pdf_out = expand(os.path.join(config["general"]["output_dir"],"pychopper/reports/normal/{{sample}}_{{unit}}_R{read}.pdf"), read=reads),
            # Reads forwarded to rescue mode
            unclass_fastq = temp(expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/normal/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)),
            # Lower-confidence reads from standard run
            rescue_fastq = temp(expand(os.path.join(config["general"]["output_dir"],"pychopper/rescued/normal/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)),
            # Reads filtered out by minimum length
            length_out = temp(expand(os.path.join(config["general"]["output_dir"], "pychopper/unclassified/normal/{{sample}}_{{unit}}_R{read}_length_out.fastq"), read=reads))
        threads: config["general"]["cores"]
        params:
            qual = config["nanopore"]["pychopqual"],
            length = config["nanopore"]["min_length"]
        conda:
            "../envs/preprocessing/pychopper.yaml"
        shell:
                """
                pychopper -m edlib -b {input.primers} -Q {params.qual} -z {params.length} -l {output.length_out} -c {input.primers_config} -r {output.pdf_out} -u {output.unclass_fastq} -w {output.rescue_fastq} -t {threads} {input.fastq} {output.out_fastq};
                """

    # Re-run pychopper with -x rescue on reads unclassified in the standard pass.
    # The rescue algorithm uses relaxed primer matching and can recover
    # reads where one primer is degraded or truncated. Reads still unclassified
    # after rescue are written to unclass_unclass_fastq for inspection.
    rule pychopper_rescue:
        input:
            # Reads unclassified in the standard pychopper run
            unclass_fastq=expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/normal/{{sample}}_{{unit}}_R{read}.fastq"), read=reads),
            # SSP/VNP primer FASTA for rescue classification
            primers = expand(os.path.join(config["general"]["output_dir"],"pychopper/custom_primers.fasta")),
            # Expected primer orientations for rescue mode
            primers_config = expand(os.path.join(config["general"]["output_dir"],"pychopper/config_primers.txt"))
        output:
            # Additional full-length reads recovered by rescue mode
            unclass_out_fastq = temp(expand(os.path.join(config["general"]["output_dir"],"pychopper/output/rescue/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)),
            # Reads still unclassified after rescue
            unclass_unclass_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/unclassified/rescue/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads),
            # Lower-confidence reads classified during rescue
            unclass_rescue_fastq = temp(expand(os.path.join(config["general"]["output_dir"],"pychopper/rescued/rescue/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)),
            # Rescue-mode classification report
            unclass_pdf = expand(os.path.join(config["general"]["output_dir"],"pychopper/reports/rescue/{{sample}}_{{unit}}_R{read}.pdf"),  read=reads),
            # Rescue reads filtered out by minimum length
            length_out = temp(expand(os.path.join(config["general"]["output_dir"], "pychopper/unclassified/rescue/{{sample}}_{{unit}}_R{read}_length_out.fastq"), read=reads))
        threads: config["general"]["cores"]
        params:
            qual =config["nanopore"]["pychopqual"],
            length = config["nanopore"]["min_length"]
        conda:
            "../envs/preprocessing/pychopper.yaml"
        shell:
            """
            pychopper -m edlib -x rescue -b {input.primers} -Q {params.qual} -z {params.length} -l {output.length_out} -c {input.primers_config} -r {output.unclass_pdf} -u {output.unclass_unclass_fastq} -w {output.unclass_rescue_fastq} -t {threads} {input.unclass_fastq} {output.unclass_out_fastq};
            """

    # Combine classified reads from the standard and rescue runs into a single
    # FASTQ file for downstream FASTA conversion and read correction.
    rule merge_pychopper:
        input:
            # Classified reads from the standard run
            out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/normal/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads),
            # Classified reads recovered by rescue mode
            unclass_out_fastq = expand(os.path.join(config["general"]["output_dir"],"pychopper/output/rescue/{{sample}}_{{unit}}_R{read}.fastq"),  read=reads)
        output:
            # Combined classified reads for downstream processing
            pychopper_merged = expand(os.path.join(config["general"]["output_dir"], "pychopper/output/merged/{{sample}}_{{unit}}_R{read}.fastq"), read=reads)
        shell:
            """
            cat {input.out_fastq} {input.unclass_out_fastq} > {output.pychopper_merged}
            """
