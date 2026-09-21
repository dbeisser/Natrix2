import os

if config['classify']['mothur']:
    rule generate_otu_fasta:
        input:
            os.path.join(config["general"]["output_dir"], "finalData/{database}/full_table.csv")
        output:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/OTU_mumu.fasta")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/generate_otu_fasta/{database}.log")
        script:
            "../scripts/utilities/generate_fasta.py"

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
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"

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
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" > {output[0]};
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

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

else:
    rule generate_otu_fasta:
        input:
            expand(os.path.join(config["general"]["output_dir"], "finalData/blast_{database}/full_table.csv"), database=config['blast']['database'].lower())
        output:
            os.path.join(config["general"]["output_dir"],"blast/OTU_mumu.fasta")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/mumu/generate_otu_fasta.log")
        script:
            "../scripts/utilities/generate_fasta.py"

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
            "vsearch --usearch_global {input} -db {input} --self --id .84 --iddef 1 " \
            "--userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"

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
            	cut -d "," -f 1,3- {input[0]} --output-delimiter="\t" > {output[0]};
            	mumu --otu_table {output[0]} --match_list {input[1]} --new_otu_table {output[1]} --log {log}
            """

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
