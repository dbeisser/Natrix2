import os
if not config['dataset']['nanopore'] and config['classify']['mothur']:
    if config["classify"]["database"] == "pr2":
        rule mothur_classify:
            input:
                os.path.join(config["general"]["output_dir"],"clustering/representatives_mod.fasta") if config["general"]["seq_rep"] == "OTU" or config['clustering'] == "swarm" else ( os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta") if config['clustering']== "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta")),
                expand("database/pr2db.{pr2_db_version}.fasta", pr2_db_version=config["database_version"]["pr2"])
            output:
                os.path.join(config["general"]["output_dir"], "mothur/pr2/mothur_out.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/pr2/mothur_out.taxonomy")
            params:
                template=config['database_path']['pr2_ref'],
                taxonomy=config['database_path']['pr2_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                threads=config['general']['cores'],
                output=config['general']['output_dir'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"]["seq_rep"] == "OTU" or config["general"]["seq_rep"] == "ASV" else os.path.join(config["general"]["output_dir"],"filtering"),
            conda:
                "../envs/mothur.yaml"
            log:
                os.path.join(config["general"]["output_dir"], "logs/mothur_classify.log")
            shell:
                """
                    mothur "#set.logfile(name={log}); classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    # sed -i "s/([^()]*)//g" {params.input}/*.taxonomy
                    input_dir=$(dirname {input[0]});
                    mv $input_dir/*.taxonomy {output[1]};
                    mv $input_dir/*.summary {output[0]};
                """

    elif config["classify"]["database"] == "unite":
        rule mothur_classify:
            input:
                os.path.join(config["general"]["output_dir"],"clustering/representatives_mod.fasta") if config["general"]["seq_rep"] == "OTU" or config['clustering'] == "swarm" else ( os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta") if config['clustering']== "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta")),
                "database/unite_v10.fasta"
            output:
                os.path.join(config["general"]["output_dir"],"mothur/unite/mothur_out.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/unite/mothur_out.taxonomy")
            params:
                template=config['database_path']['unite_ref'],
                taxonomy=config['database_path']['unite_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                output=config['general']['output_dir'],
                threads=config['general']['cores'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"]["seq_rep"] == "OTU"   or config["general"]["seq_rep"] == "ASV" else os.path.join(config["general"]["output_dir"],"filtering"),
            conda:
                "../envs/mothur.yaml"
            log:
                "logs/mothur_classify.log"
            shell:
                """
                    mothur "#classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    # sed -i "s/([^()]*)//g" {params.input}/*.taxonomy
                    input_dir=$(dirname {input[0]});
                    mv $input_dir/*.taxonomy {output[1]};
                    mv $input_dir/*.summary {output[0]};
                """

    elif config["classify"]["database"] == "silva":
        rule mothur_classify:
            input:
                os.path.join(config["general"]["output_dir"],"clustering/representatives_mod.fasta") if config["general"]["seq_rep"] == "OTU" or config['clustering'] == "swarm" else ( os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta") if config['clustering']== "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta")),
                expand(["database/silva_db.{silva_db_version}.fasta", "database/silva_db.{silva_db_version}.tax"], silva_db_version=config["database_version"]["silva"])
            output:
                os.path.join(config["general"]["output_dir"],"mothur/silva/mothur_out_raw.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/silva/mothur_out_raw.taxonomy")
            params:
                template=config['database_path']['silva_ref'],
                taxonomy=config['database_path']['silva_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                output=config['general']['output_dir'],
                threads=config['general']['cores'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"]["seq_rep"] == "OTU"  or config["general"]["seq_rep"] == "ASV" else os.path.join(config["general"]["output_dir"],"filtering"),
            conda:
                "../envs/mothur.yaml"
            log:
                os.path.join(config["general"]["output_dir"],"logs/mothur_classify.log")
            shell:
                """
                    mothur "#classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    # sed -i "s/([^()]*)//g" {params.input}/*.taxonomy
                    input_dir=$(dirname {input[0]});
                    mv $input_dir/*.taxonomy {output[1]};
                    mv $input_dir/*.summary {output[0]};
                """
        
        # Removes redundant 'unclassified' entries from the taxonomy file (/mothur_out_raw.taxonomy)
        rule filter_unclassified:
            input:
                raw_taxonomy=os.path.join(config["general"]["output_dir"], "mothur/silva/mothur_out_raw.taxonomy"),
                raw_summary=os.path.join(config["general"]["output_dir"], "mothur/silva/mothur_out_raw.summary")
            output:
                cleaned_taxonomy=os.path.join(config["general"]["output_dir"], "mothur/silva/mothur_out.taxonomy"),
                cleaned_summary=os.path.join(config["general"]["output_dir"], "mothur/silva/mothur_out.summary")
            script:
                "../scripts/filter_unclassified.py"

    ##
    rule merge_output:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/mothur_out.taxonomy"),
            os.path.join(config["general"]["output_dir"], "clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and config['clustering']=="swarm" else (os.path.join(config["general"]["output_dir"], "clustering/vsearch_table.csv") if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"], "filtering/filtered_table.csv")),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/metadata_table.csv")
        params:
            clustering=config['clustering']
        script:
                "../scripts/merge_results2.py"

######################################################################################################################
elif config['dataset']['nanopore']:
    if config["classify"]["database"] == "pr2":
        rule mothur_classify:
            input:
                expand(os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta")) if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
                expand("database/pr2db.{pr2_db_version}.fasta",pr2_db_version=config["database_version"]["pr2"])
            output:
                os.path.join(config["general"]["output_dir"],"mothur/pr2/mothur_out.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/pr2/mothur_out.taxonomy")
            params:
                template=config['database_path']['pr2_ref'],
                taxonomy=config['database_path']['pr2_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                threads=config['general']['cores'],
                output=config['general']['output_dir'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"][
                                                                                        "seq_rep"] == "OTU" and config ["dataset"]["nanopore"] == "FALSE" 
                   else ( os.path.join(config["general"]["output_dir"],"clustering") if config["clustering"]=="vsearch" else os.path.join(config["general"]["output_dir"],"filtering")),
            conda:
                "../envs/mothur.yaml"
            log:
                os.path.join(config["general"]["output_dir"],"logs/mothur_classify.log")
            shell:
                """
                    mothur "#set.logfile(name={log}); classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    #sed -i "s/([^()]*)//g" {params.input}/*.taxonomy 
                    mv {params.input}/*.taxonomy {params.output}/mothur/pr2/mothur_out.taxonomy;
                    mv {params.input}/*.summary {params.output}/mothur/pr2/mothur_out.summary;
                """

    elif config["classify"]["database"] == "unite":
        rule mothur_classify:
            input:
                expand(os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta")) if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
                "database/unite_v10.fasta"
            output:
                os.path.join(config["general"]["output_dir"],"mothur/unite/mothur_out.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/unite/mothur_out.taxonomy")
            params:
                template=config['database_path']['unite_ref'],
                taxonomy=config['database_path']['unite_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                output=config['general']['output_dir'],
                threads=config['general']['cores'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"][
                                                                                        "seq_rep"] == "OTU" and config ["dataset"]["nanopore"] == "FALSE" 
                   else ( os.path.join(config["general"]["output_dir"],"clustering") if config["clustering"]=="vsearch" else os.path.join(config["general"]["output_dir"],"filtering")),
            conda:
                "../envs/mothur.yaml"
            log:
                os.path.join(config["general"]["output_dir"],"logs/mothur_classify.log")
            shell:
                """
                    mothur "#classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    #sed -i "s/([^()]*)//g" {params.input}/*.taxonomy 
                    mv {params.input}/*.taxonomy {params.output}/mothur/unite/mothur_out.taxonomy;
                    mv {params.input}/*.summary {params.output}/mothur/unite/mothur_out.summary;
                """

    elif config["classify"]["database"] == "silva":
        rule mothur_classify:
            input:
                expand(os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta")) if config['clustering']== "vsearch" else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
                expand(["database/silva_db.{silva_db_version}.fasta",
                "database/silva_db.{silva_db_version}.tax"],silva_db_version=config["database_version"]["silva"])
            output:
                os.path.join(config["general"]["output_dir"],"mothur/silva/mothur_out.summary"),
                os.path.join(config["general"]["output_dir"],"mothur/silva/mothur_out.taxonomy")
            params:
                template=config['database_path']['silva_ref'],
                taxonomy=config['database_path']['silva_tax'],
                search=config['classify']['search'],
                method=config['classify']["method"],
                cutoff=config['classify']['cutoff'],
                output=config['general']['output_dir'],
                threads=config['general']['cores'],
                input=os.path.join(config["general"]["output_dir"],"clustering") if config["general"][
                                                                                        "seq_rep"] == "OTU" and config ["dataset"]["nanopore"] == "FALSE" 
                   else ( os.path.join(config["general"]["output_dir"],"clustering") if config["clustering"]=="vsearch" else os.path.join(config["general"]["output_dir"],"filtering")),
            conda:
                "../envs/mothur.yaml"
            log:
                os.path.join(config["general"]["output_dir"],"logs/mothur_classify.log")
            shell:
                """
                    mothur "#classify.seqs(fasta={input[0]}, cutoff={params.cutoff}, reference={params.template}, taxonomy={params.taxonomy}, method={params.method}, processors={params.threads}, output=simple, search={params.search})";
                    #sed -i "s/([^()]*)//g" {params.input}/*.taxonomy 
                    mv {params.input}/*.taxonomy {params.output}/mothur/silva/mothur_out.taxonomy;
                    mv {params.input}/*.summary {params.output}/mothur/silva/mothur_out.summary;
                """

    rule merge_output:
        input:
            os.path.join(config["general"]["output_dir"],"mothur/{database}/mothur_out.taxonomy"),
            expand(os.path.join(config["general"]["output_dir"], "clustering/vsearch_table.csv")) if config['clustering'] == "vsearch" else os.path.join(config["general"]["output_dir"], "filtering/filtered_table.csv"),
        output:
            os.path.join(config["general"]["output_dir"],"finalData/{database}/full_table.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/OTU_table.csv"),
            os.path.join(config["general"]["output_dir"],"finalData/{database}/metadata_table.csv")
        params:
            clustering=config['clustering']
        script:
            "../scripts/merge_results2.py"


