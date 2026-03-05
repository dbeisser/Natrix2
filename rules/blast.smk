import os

if config["blast"]["database"] == "SILVA":

    rule make_silva_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"]),
            config["blast"]["db_path"] + ".fasta"
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            """
                dir_name=$(dirname {params[0]});
                latest_file=$(wget -qO- https://ftp.arb-silva.de/current/Exports/ | grep -oP 'SILVA_\d+\.\d+_SSURef_tax_silva.fasta.gz' | sort -V | tail -n 1);
                wget -P $dir_name/ --progress=bar https://ftp.arb-silva.de/current/Exports/$latest_file;
                pigz -d -c $dir_name/$latest_file > $dir_name/silva.db.fasta;
                makeblastdb -in $dir_name/silva.db.fasta -dbtype nucl -parse_seqids -out $dir_name/silva.db -blastdb_version 5
            """

    rule create_silva_taxonomy:
        input: config["blast"]["db_path"] + ".fasta"
        output: os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        conda:
            "../envs/blast.yaml"
        script:
            "../scripts/create_silva_taxonomy.py"

elif config["blast"]["database"] == "NCBI":

    rule make_ncbi_db:
        output:
            marker=config["blast"]["db_path"] + config["blast"]["db_type"] + ".downloaded.ok"
        params:
            db_dir=config["blast"]["db_path"],
            db_type=config["blast"]["db_type"]
        log:
            config["blast"]["db_path"] + config["blast"]["db_type"] + ".download.log"
        conda:
            "../envs/blast.yaml"
        shell:
            """
                mkdir -p {params.db_dir}
                echo "[{params.db_type}] Downloading into: {params.db_dir}" > {log}
                ( cd {params.db_dir} && update_blastdb.pl --source aws --decompress --verbose --verbose {params.db_type} ) >> {log} 2>&1
                echo "[{params.db_type}] Done." >> {log}
                touch {output.marker}
                sleep 10
            """

    rule download_taxonomy:
        output:
            os.path.join(config["blast"]["db_path"], "fullnamelineage.dmp")
        params:
            db_path=config["blast"]["db_path"]
        conda:
            "../envs/blast.yaml"
        shell:
            """
                dir_name={params.db_path}
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
                tar xvzf $dir_name/new_taxdump.tar.gz -C $dir_name
            """

    rule create_blast_taxonomy:
        input:
            os.path.join(config["blast"]["db_path"], "fullnamelineage.dmp")
        output:
            os.path.join(config["blast"]["db_path"], "tax_lineage.h5")
        conda:
            "../envs/blast.yaml"
        log:
            os.path.join(config["general"]["output_dir"],"logs/BLAST.log")
        script:
            "../scripts/create_blast_taxonomy.py"


rule blast:
    input:
        query = os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta") \
            if config["general"]["seq_rep"] == "OTU" and not config['dataset']['nanopore'] \
            else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),

        db_ready = (
            config["blast"]["db_path"] + config["blast"]["db_type"] + ".downloaded.ok"
            if config["blast"]["database"] == "NCBI"
            else expand(config["blast"]["db_path"] + "{ext}",
                        ext=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"])
            if config["blast"]["database"] == "SILVA"
            else []
        )
    output:
        os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv")
    log:
        os.path.join(config["general"]["output_dir"],"blast/blast.log")
    threads: 
        config["general"]["cores"]
    params:
        db_path = (
            config["blast"]["db_path"] + config["blast"]["db_type"]
            if config["blast"]["database"] == "NCBI"
            else config["blast"]["db_path"]
        ),
        max_target_seqs=str(config["blast"]["max_target_seqs"]) if config["blast"]["database"] == "NCBI" else "1",
        ident=str(config["blast"]["ident"]),
        evalue=str(config["blast"]["evalue"]),
        out6='"6 qseqid qlen length pident mismatch qstart qend sstart send gaps evalue staxid sseqid"'
    conda:
        "../envs/blast.yaml"
    shell:
        """
            echo "Running BLAST with database: {params.db_path}" > {log}
            echo "Query file: {input.query}" >> {log}
            echo "Threads: {threads}" >> {log}
            echo "Database info:" >> {log}
            blastdbcmd -db {params.db_path} -info >> {log} 2>&1

            blastn -num_threads {threads} \
                -query {input.query} \
                -db {params.db_path} \
                -max_target_seqs {params.max_target_seqs} \
                -perc_identity {params.ident} \
                -evalue {params.evalue} \
                -outfmt {params.out6} \
                -out {output} \
                2>> {log}
        """

if config['blast']['blast']:
	if config["blast"]["database"]=="NCBI":

		rule ncbi_taxonomy:
			input:
				blast_result = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv"),
                lineage = os.path.join(config["blast"]["db_path"], "tax_lineage.h5")
			output:
				tax_lineage = temp(os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")),
				all_tax = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy_all.tsv")
			params:
				max_target_seqs=config["blast"]["max_target_seqs"],
				drop_tax_classes=str(config["blast"]["drop_tax_classes"])
			conda:
				"../envs/blast.yaml"
			log:
				os.path.join(config["general"]["output_dir"],"logs/BLAST.log")
			script:
				"../scripts/ncbi_taxonomy.py"

		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and not config['dataset']['nanopore'] else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			conda:
				"../envs/merge_results.yaml"
			log:
				os.path.join(config["general"]["output_dir"],"logs/BLAST.log")
			script:
				"../scripts/merge_results.py"


	elif config["blast"]["database"] == "SILVA":

		rule silva_taxonomy:
			input:
				blast_result = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv"),
				lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
			output:
				temp(os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv"))
			params:
				drop_tax_classes=str(config["blast"]["drop_tax_classes"])
			conda:
				"../envs/blast.yaml"
			log:
				os.path.join(config["general"]["output_dir"],"logs/BLAST.log")
			script:
				"../scripts/silva_taxonomy.py"

		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") if config["general"]["seq_rep"] == "OTU" and not config['dataset']['nanopore']  else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			conda:
				"../envs/merge_results.yaml"
			log:
				os.path.join(config["general"]["output_dir"],"logs/BLAST.log")
			script:
				"../scripts/merge_results.py"
