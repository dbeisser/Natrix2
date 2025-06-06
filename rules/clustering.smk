import os
if not config['dataset']['nanopore'] and config['clustering']== 'swarm':
    # dada runs per sample after preprocessing and before chimera removal and AmpliconDuo
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
        conda:
            "../envs/dada2.yaml"
        log:
            os.path.join(config["general"]["output_dir"], "logs/{sample}_{unit}_dada2.log")
        script:
            "../scripts/dada2.R"


    rule swarm:
        input:
            os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta"),
            os.path.join(config["general"]["output_dir"],"clustering/merged.swarms")
        threads: config["general"]["cores"]
        conda:
            "../envs/swarm.yaml"
        shell:
            "swarm -t {threads} -f -z -w {output[0]} < {input} > {output[1]}"

    rule swarm_headers:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta"),
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"clustering/representatives_mod.fasta")
        script:
            "../scripts/swarm_fasta_headers.py"


    rule clust_merge_results:
        input:
            merged=os.path.join(config["general"]["output_dir"],"clustering/merged.swarms") if config['clustering'] == "swarm" else os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt"),
            final_table_path2=os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
        output:
            out=os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config['clustering'] == "swarm" else os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
        conda:
            "../envs/merge_results.yaml"
        script:
            "../scripts/merge_clust_results.py"

    # SWARM clustering runs on all samples after AmpliconDuo
    rule write_fasta:
        input:
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv")
        output:
            os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
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
