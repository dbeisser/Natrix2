import os

# Prepare per-sample, per-unit read files using the configured demultiplexing script.
# The generated compressed FASTQ files are temporary workflow intermediates.
rule demultiplex:
    output:
        temp(expand(os.path.join(config["general"]["output_dir"],"demultiplexed/{unit.sample}_{unit.unit}_R{read}.fastq.gz"), unit=units.reset_index().itertuples(), read=reads))
    params:
        filename = config["general"]["filename"].rstrip("/"),
        primertable = config["general"]["primertable"],
        # demultiplexing = config['general']['demultiplexing'],
        read_sorting = config['general']['read_sorting'],
        assembled = config['general']['already_assembled'],
        name_ext = config['merge']['name_ext'],
        output_dir = config['general']['output_dir']
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/demultiplexing/demultiplex.log")
    conda:
        "../envs/preprocessing/demultiplexing.yaml"
    script:
        "../scripts/preprocessing/demultiplexing.py"

# Decompress each demultiplexed FASTQ file into a temporary uncompressed file.
rule unzip:
    input:
         os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_R{read}.fastq.gz")
    output:
        temp(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.tmp"))
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/demultiplexing/unzip/{sample}_{unit}_{read}.log")
    shell:
        # Decompress the input with pigz and write the FASTQ data to the output file.
        # Only stderr is logged (file and terminal); the decompressed data go to the output file.
         "pigz -d -c {input} 2>&1 > {output} | tee {log:q}"

# Convert CRLF line endings to LF when detected; otherwise rename the FASTQ file.
rule check_format:
    input:
        os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.tmp")
    output:
        temp(os.path.join(config["general"]["output_dir"],"demultiplexed/{sample}_{unit}_{read}.fastq"))
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/demultiplexing/check_format/{sample}_{unit}_{read}.log")
    shell:
        """
            # Inspect the input file and check whether its line endings contain CRLF.
            if find {input} -not -type d -exec file '{{}}' ';' | grep CRLF
            then
                sed 's/\r$//' {input} 2>&1 > {output} | tee {log:q}
            else
                mv {input} {output} 2>&1 | tee {log:q}
            fi
        """
