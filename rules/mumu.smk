import os

# Curate OTU tables using VSEARCH sequence similarities and MUMU.
# Prepare sequence matches and abundance tables, then combine the curated
# results with existing classification tables.

# [Classifier-specific workflow: MOTHUR]
if config['classify']['mothur']:
    # Extract sequences from the MOTHUR result table to a FASTA file.
    rule generate_otu_fasta:
        input:
            os.path.join(config["general"]["output_dir"], "finalData/{database}/full_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_mumu.fasta")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/generate_otu_fasta/{database}.log")
        script:
            "../scripts/utilities/generate_fasta.py"

    # Search OTU sequences against the same FASTA using VSEARCH.
    # Require at least 84% identity and 90% query coverage; write query,
    # target, and identity fields for downstream MUMU processing.
    rule vsearch_otu:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_mumu.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/match_scores.txt")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/vsearch_otu/{database}.log")
        conda:
            "../envs/analysis/vsearch.yaml"
        shell:
            # Run VSEARCH and mirror diagnostic output to the console and log.
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10" \
            " 2>&1 | tee {log:q}"

    # Convert the abundance CSV to a tab-separated MUMU input table,
    # excluding its second column; run MUMU using the VSEARCH match list.
    rule run_mumu:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config ['clustering']== "swarm" or config["dataset"]["nanopore"] == "FALSE"  else (os.path.join(config["general"]["output_dir"], "clustering/vsearch_table.csv") if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"], "filtering/filtered_table.csv")),
            expand(os.path.join(config["general"]["output_dir"],"mothur/{database}/match_scores.txt"), database=config['classify']['database']),
        output:
            temp(expand(os.path.join(config["general"]["output_dir"], "mothur/{database}/OTU_table_mumu.tmp"), database=config['classify']['database'])),
            expand(os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_table_mumu.csv"), database=config['classify']['database'])
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/run_mumu.log")
        conda:
            "../envs/classification/mumu.yaml"
        shell:
            """
            	# Write the MUMU input table; send diagnostic output through tee.
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" 2>&1 > {output[0]} | tee {log:q};
            	# Run MUMU and request its tool log using --log.
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

    # Combine the curated MUMU table with the existing MOTHUR result table
    # and write the resulting data, abundance, and metadata tables.
    rule merge_mumu_mothur_output:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_table_mumu.csv"),
        	os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table.csv"),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/metadata_table_mumu.csv")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/merge_mumu_mothur_output/{database}.log")
        script:
            "../scripts/utilities/merge_mumu_output.py"

# [Classifier-specific workflow: BLAST]
else:
    # Extract sequences from the BLAST result table to a FASTA file.
    rule generate_otu_fasta:
        input:
            expand(os.path.join(config["general"]["output_dir"], "finalData/blast_{database}/full_table.csv"), database=config['blast']['database'].lower())
        output:
            os.path.join(config["general"]["output_dir"],"blast/OTU_mumu.fasta")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/generate_otu_fasta.log")
        script:
            "../scripts/utilities/generate_fasta.py"

    # Search OTU sequences against the same FASTA using VSEARCH.
    # Require at least 84% identity and 90% query coverage; write query,
    # target, and identity fields for downstream MUMU processing.
    rule vsearch_otu:
        input:
            os.path.join(config["general"]["output_dir"],"blast/OTU_mumu.fasta")
        output:
            os.path.join(config["general"]["output_dir"],"blast/match_scores.txt")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/vsearch_otu.log")
        conda:
            "../envs/analysis/vsearch.yaml"
        shell:
            # Run VSEARCH and mirror diagnostic output to the console and log.
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10" \
            " 2>&1 | tee {log:q}"

    # Convert the abundance CSV to a tab-separated MUMU input table,
    # excluding its second column; run MUMU using the VSEARCH match list.
    rule run_mumu:
        input:
            os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config ['clustering']== "swarm" or config["dataset"]["nanopore"] == "FALSE"  else (os.path.join(config["general"]["output_dir"], "clustering/vsearch_table.csv") if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"], "filtering/filtered_table.csv")),
            expand(os.path.join(config["general"]["output_dir"],"blast/match_scores.txt"), database=config['classify']['database']),
        output:
            temp(expand(os.path.join(config["general"]["output_dir"], "blast/OTU_table_mumu.tmp"), database=config['classify']['database'])),
            expand(os.path.join(config["general"]["output_dir"],"blast/OTU_table_mumu.csv"), database=config['classify']['database'])
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/run_mumu.log")
        conda:
            "../envs/classification/mumu.yaml"
        shell:
            """
            	# Write the MUMU input table; send diagnostic output through tee.
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" 2>&1 > {output[0]} | tee {log:q};
            	# Run MUMU and request its tool log using --log.
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

    # Combine the curated MUMU table with the existing BLAST result table.
    # Write result tables and a file for sequence IDs not merged by the script.
    rule merge_mumu_blast_output:
        input:
            os.path.join(config["general"]["output_dir"],"blast/OTU_table_mumu.csv"),
        	os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/full_table.csv"),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/full_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/OTU_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/metadata_table_mumu.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/blast_{database}/unmerged_seqids.csv")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/merge_mumu_blast_output/{database}.log")
        script:
            "../scripts/utilities/merge_mumu_output_2.py"
