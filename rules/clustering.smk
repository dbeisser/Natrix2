import os

# Process Illumina reads and prepare sequence clusters for classification.
# This section defines DADA2 processing, FASTA preparation, SWARM clustering,
# cluster/abundance-table merging, and representative FASTA header updates.
# These rules are defined only for the non-Nanopore SWARM configuration.
if not config['dataset']['nanopore'] and config['clustering']== 'swarm':

    # Process per-sample trimmed reads using the DADA2 R script.
    # Pass paired-end settings, minimum overlap, and filtering configuration.
    # Protect the resulting FASTA against accidental removal by Snakemake.
    rule DADA2:
        input:
            fwd=os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_1_cut.fastq"),
            rev=os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_2_cut.fastq") if config["merge"]["paired_End"] == True else []
        output:
            protected(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_dada.fasta"))
        params:
            paired_end=config["merge"]["paired_End"],
            minoverlap=config["qc"]["minoverlap"],
            splitsamples=config["merge"]["filter_method"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/clustering/DADA2/{sample}_{unit}.log")
        conda:
            "../envs/classification/dada2.yaml"
        script:
            "../scripts/classification/dada2.R"

    # Cluster filtered FASTA sequences with SWARM in fastidious mode (-f).
    # Use abundance annotations (-z), write representatives (-w), and save
    # cluster membership from stdout; capture tool diagnostics in the log.
    rule swarm:
        input:
            os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta"),
            os.path.join(config["general"]["output_dir"],"clustering/merged.swarms")
        threads: config["general"]["cores"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/clustering/swarm.log")
        conda:
            "../envs/analysis/swarm.yaml"
        shell:
            # Run SWARM: cluster list (stdout) goes to merged.swarms, stderr is logged to file and terminal
            "swarm -t {threads} -f -z -w {output[0]} < {input} 2>&1 > {output[1]} | tee {log:q}"

    # Update representative FASTA headers using the SWARM abundance table.
    rule swarm_headers:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta"),
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"clustering/representatives_mod.fasta")
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/clustering/swarm_headers.log")
        script:
            "../scripts/analysis/swarm_fasta_headers.py"

    # Combine SWARM cluster memberships with the filtered abundance table.
    rule clust_merge_results:
        input:
            merged=os.path.join(config["general"]["output_dir"],"clustering/merged.swarms") if config['clustering'] == "swarm" else os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt"),
            final_table_path2=os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
        output:
            out=os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config['clustering'] == "swarm" else os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/clustering/clust_merge_results.log")
        conda:
            "../envs/utilities/merge_results.yaml"
        script:
            "../scripts/utilities/merge_clust_results.py"

    # Write filtered sequences to FASTA and add total abundances to their IDs.
    # Input rows are converted into numbered sequence identifiers.
    rule write_fasta:
        input:
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv")
        output:
            os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/clustering/write_fasta.log")
        run:
            import csv
            with open(input[0],"r") as csv_in, open(
                    output[0],"w") as fasta_out, open(output[1],"w") as csv_out:
                filtered_table = csv.reader(csv_in)
                filtered_table_seqid = csv.writer(csv_out)
                for row in enumerate(filtered_table):
                    if row[0] != 0:
                        abu_sum = sum([int(num) for num in row[1][1:]])
                        filtered_table_seqid.writerow([">{};size={};".format(
                            row[0],abu_sum)] + row[1])
                        fasta_out.write(">{};size={};\n{}\n".format(row[0],
                            abu_sum,row[1][0]))
                    else:
                        filtered_table_seqid.writerow(["seqid"] + row[1])
