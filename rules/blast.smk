import os

# Assign taxonomy to query sequences using BLAST and reference lineage data.
# Depending on configuration, the workflow prepares SILVA, NCBI, ROD, or
# Eukaryome BLAST database resources and a taxonomy lookup table.
# BLASTN searches the selected query FASTA against the database.
# Subsequent rules annotate hits with lineage information and merge
# annotations with an abundance table for downstream analysis.

# [Database setup: SILVA]
if config["blast"]["database"] == "silva":

    # Retrieve a SILVA SSU Ref FASTA release and build a nucleotide BLAST database.
    # Select the highest versioned matching file listed in the SILVA FTP directory.
    # Decompress the FASTA with pigz and index it using makeblastdb
    # with BLAST database format version 5.
    rule make_silva_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"]),
            config["blast"]["db_path"] + ".fasta"
        params:
            db_path=config["blast"]["db_path"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/make_silva_db.log")
        conda:
            "../envs/classification/blast.yaml"
        shell:
            """
                # Find the latest SILVA release (variables are not piped, otherwise they would be lost)
                dir_name=$(dirname {params[0]})
                latest_file=$(wget -qO- https://ftp.arb-silva.de/current/Exports/ | grep -oP 'SILVA_\d+\.\d+_SSURef_tax_silva.fasta.gz' | sort -V | tail -n 1)
                echo "Latest SILVA release: $latest_file" 2>&1 | tee {log:q}

                # Download and decompress SILVA (stderr to log, sequences to FASTA file)
                wget -P $dir_name/ --progress=bar https://ftp.arb-silva.de/current/Exports/$latest_file 2>&1 | tee -a {log:q}
                pigz -d -c $dir_name/$latest_file 2>&1 > $dir_name/silva.db.fasta | tee -a {log:q}

                # Build BLAST database
                makeblastdb -in $dir_name/silva.db.fasta -dbtype nucl -parse_seqids -out $dir_name/silva.db -blastdb_version 5 2>&1 | tee -a {log:q}
            """

    # Build a taxonomic lineage table from the SILVA reference FASTA headers.
    # The lineage information is stored as an HDF5 file (tax_lineage.h5),
    # which is later used to translate BLAST hits into full taxonomic ranks.
    rule create_silva_taxonomy:
        input: 
            config["blast"]["db_path"] + ".fasta"
        output: 
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/create_silva_taxonomy.log")
        conda:
            "../envs/classification/blast.yaml"
        script:
            "../scripts/classification/create_silva_taxonomy.py"


# [Database setup: NCBI]
elif config["blast"]["database"] == "ncbi":

    # Download the configured preformatted NCBI BLAST database from AWS
    # using update_blastdb.pl with decompression enabled.
    # Record completion using a marker file
    # (.downloaded.ok), which is the declared rule output.
    rule make_ncbi_db:
        output:
            marker=config["blast"]["db_path"] + config["blast"]["db_type"] + ".downloaded.ok"
        params:
            db_dir=config["blast"]["db_path"],
            db_type=config["blast"]["db_type"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/make_ncbi_db.log")
        conda:
            "../envs/classification/blast.yaml"
        shell:
            """
                # Download NCBI database; log stdout and stderr to file and terminal
                mkdir -p {params.db_dir} 2>&1 | tee {log:q}
                echo "[{params.db_type}] Downloading into: {params.db_dir}" 2>&1 | tee -a {log:q}
                ( cd {params.db_dir} && update_blastdb.pl --source aws --decompress --verbose --verbose {params.db_type} ) 2>&1 | tee -a {log:q}
                echo "[{params.db_type}] Done." 2>&1 | tee -a {log:q}

                # Create marker file after successful download
                touch {output.marker}
                sleep 10
            """

    # Download the NCBI taxonomy dump (new_taxdump) and the BLAST taxonomy
    # database (taxdb). The taxonomy dump is extracted to provide
    # fullnamelineage.dmp, which links NCBI taxonomy IDs to full lineages.
    rule download_taxonomy:
        output:
            os.path.join(config["blast"]["db_path"], "fullnamelineage.dmp")
        params:
            db_path=config["blast"]["db_path"]
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/download_taxonomy.log")
        conda:
            "../envs/classification/blast.yaml"
        shell:
            """
                # Download NCBI taxonomy; log stdout and stderr to file and terminal
                dir_name={params.db_path}
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz 2>&1 | tee {log:q}
                wget -N -P $dir_name/ --progress=bar ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz 2>&1 | tee -a {log:q}
                tar xvzf $dir_name/new_taxdump.tar.gz -C $dir_name 2>&1 | tee -a {log:q}
            """

    # Convert the NCBI lineage file (fullnamelineage.dmp) into an HDF5 lookup
    # table (tax_lineage.h5) for fast mapping of BLAST taxonomy IDs to lineages.
    rule create_blast_taxonomy:
        input:
            os.path.join(config["blast"]["db_path"], "fullnamelineage.dmp")
        output:
            os.path.join(config["blast"]["db_path"], "tax_lineage.h5")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/create_blast_taxonomy.log")
        conda:
            "../envs/classification/blast.yaml"
        script:
            "../scripts/classification/create_blast_taxonomy.py"


# [Database setup: ROD]
elif config["blast"]["database"] == "rod":

    # Download the ROD reference database release from GitHub and build a
    # BLAST database. The FASTA headers are cleaned ('|' replaced by spaces),
    # indexed with makeblastdb, and temporary download files are removed.
    rule make_ROD_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"]),
            config["blast"]["db_path"] + ".fasta"
        params:
            db_path=config["blast"]["db_path"],
            db_version=config["database_version"]["rod"],
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/make_ROD_db.log")
        conda:
            "../envs/classification/blast.yaml"
        shell:
            """
                # Download and extract ROD; log stdout and stderr to file and terminal
                wget -O database/v{params.db_version}-goldenrod.zip --progress=bar https://github.com/krabberod/ROD/archive/refs/tags/v{params.db_version}-goldenrod.zip 2>&1 | tee {log:q}
                unzip -o database/v{params.db_version}-goldenrod.zip -d database/ 2>&1 | tee -a {log:q}

                # Decompress and clean FASTA headers (stderr to log, sequences to FASTA file)
                gunzip -c database/ROD-{params.db_version}-goldenrod/ROD_v{params.db_version}_reference_sequences.fasta.gz 2>&1 > {params.db_path}.{params.db_version}_reference_sequences_lang.fasta | tee -a {log:q}
                sed 's/|/ /g' {params.db_path}.{params.db_version}_reference_sequences_lang.fasta 2>&1 > {params.db_path}.{params.db_version}_reference_sequences.fasta | tee -a {log:q}

                # Build BLAST database
                makeblastdb -in {params.db_path}.{params.db_version}_reference_sequences.fasta -dbtype nucl -parse_seqids -out {params.db_path} -blastdb_version 5 2>&1 | tee -a {log:q}

                # Rename FASTA and remove temporary files
                mv {params.db_path}.{params.db_version}_reference_sequences.fasta {params.db_path}.fasta 2>&1 | tee -a {log:q}
                rm database/v{params.db_version}-goldenrod.zip 2>&1 | tee -a {log:q}
                rm -r database/ROD-{params.db_version}-goldenrod 2>&1 | tee -a {log:q}
            """

    # Build a taxonomic lineage table from the ROD reference FASTA headers
    # and store it as an HDF5 file (tax_lineage.h5) for later BLAST annotation.
    rule create_ROD_taxonomy:
        input: 
            config["blast"]["db_path"] + ".fasta"
        output: 
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/create_ROD_taxonomy.log")
        conda:
            "../envs/classification/blast.yaml"
        script:
            "../scripts/classification/create_rod_taxonomy.py"


# [Database setup: Eukaryome]
elif config["blast"]["database"] == "eukaryome":

    # Download the Eukaryome reference database (mothur format) and build a
    # BLAST database. The FASTA and taxonomy (.tax) files are extracted,
    # the FASTA is indexed with makeblastdb, and temporary files are removed.
    rule make_Eukaryome_db:
        output:
            expand(config["blast"]["db_path"] + "{file_extension}", file_extension=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"]),
            config["blast"]["db_path"] + ".fasta",
            config["blast"]["db_path"] + ".tax"
        params:
            db_path=config["blast"]["db_path"],
            eukaryome_version=config["database_path"]["eukaryome_version"],
            db_version=config["database_version"]["eukaryome"],
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/make_Eukaryome_db.log")
        conda:
            "../envs/classification/blast.yaml"
        shell:
            """
                # Download and extract Eukaryome; log stdout and stderr to file and terminal
                wget -P ./database/ --progress=bar https://sisu.ut.ee/wp-content/uploads/sites/643/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip 2>&1 | tee {log:q}
                python3 -m zipfile -e database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip database/ 2>&1 | tee -a {log:q}

                # Copy FASTA and taxonomy files (stderr to log, content to file)
                cat database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.fasta 2>&1 > {params.db_path}.{params.db_version}.fasta | tee -a {log:q}
                cat database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.tax 2>&1 > {params.db_path}.tax | tee -a {log:q}

                # Build BLAST database
                makeblastdb -in {params.db_path}.{params.db_version}.fasta -dbtype nucl -parse_seqids -out {params.db_path} -blastdb_version 5 2>&1 | tee -a {log:q}

                # Remove temporary files and rename FASTA
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip 2>&1 | tee -a {log:q}
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.fasta 2>&1 | tee -a {log:q}
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.tax 2>&1 | tee -a {log:q}
                mv {params.db_path}.{params.db_version}.fasta {params.db_path}.fasta 2>&1 | tee -a {log:q}
            """

    # Build a taxonomic lineage table from the Eukaryome taxonomy file (.tax)
    # and store it as an HDF5 file (tax_lineage.h5) for later BLAST annotation.
    rule create_Eukaryome_taxonomy:
        input: 
            config["blast"]["db_path"] + ".tax"
        output: 
            os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
        log:
            os.path.join(config["general"]["output_dir"],"logfiles/blast/create_Eukaryome_taxonomy.log")
        conda:
            "../envs/classification/blast.yaml"
        script:
            "../scripts/classification/create_eukaryome_taxonomy.py"


# Align the representative sequences against the selected BLAST database.
# The query file depends on the dataset and clustering method (OTU
# representatives, vsearch centroids, or filtered sequences).
# BLASTN uses the configured identity and E-value thresholds. For NCBI,
# the configured maximum hit count is used; other databases use a limit of one.
# Results are written as a tabular BLAST output (outfmt 6).
rule blast:
    input:
        # Query sequences: OTU representatives (Illumina), vsearch centroids
        # (Nanopore + vsearch) or filtered sequences (all other cases)
        query = os.path.join(config["general"]["output_dir"],"clustering/representatives.fasta") \
            if config["general"]["seq_rep"] == "OTU" and not config["dataset"]["nanopore"] \
            else os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta") \
            if config["clustering"]== "vsearch" and config["dataset"]["nanopore"] \
            else os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),

        # Database files: marker file for NCBI, index files for local databases
        db_ready = (config["blast"]["db_path"] + config["blast"]["db_type"] + ".downloaded.ok" if config["blast"]["database"] == "ncbi"
            else expand(config["blast"]["db_path"] + "{ext}", ext=[".ndb", ".nhr", ".nin", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto"])
            if config["blast"]["database"] in ["silva", "rod", "eukaryome"] else [])
    output:
        os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv")
    params:
        db_path = (config["blast"]["db_path"] + config["blast"]["db_type"]
            if config["blast"]["database"] == "ncbi"
            else config["blast"]["db_path"]),
        max_target_seqs=str(config["blast"]["max_target_seqs"]) if config["blast"]["database"] == "ncbi" else "1",
        ident=str(config["blast"]["ident"]),
        evalue=str(config["blast"]["evalue"]),
        out6='"6 qseqid qlen length pident mismatch qstart qend sstart send gaps evalue staxid sseqid"'
    threads: 
        config["general"]["cores"]
    log:
        os.path.join(config["general"]["output_dir"],"logfiles/blast/blast.log")
    conda:
        "../envs/classification/blast.yaml"
    shell:
        """
            # Write run information to log
            echo "Running BLAST with database: {params.db_path}" 2>&1 | tee {log:q}
            echo "Query file: {input.query}" 2>&1 | tee -a {log:q}
            echo "Threads: {threads}" 2>&1 | tee -a {log:q}
            echo "Database info:" 2>&1 | tee -a {log:q}
            blastdbcmd -db {params.db_path} -info 2>&1 | tee -a {log:q}

            # Run BLAST
            blastn \
                -num_threads {threads} \
                -query {input.query} \
                -db {params.db_path} \
                -max_target_seqs {params.max_target_seqs} \
                -perc_identity {params.ident} \
                -evalue {params.evalue} \
                -outfmt {params.out6} \
                -out {output} \
                2>&1 | tee -a {log:q}
        """


# [Taxonomy assignment and result merging for the selected database]
if config['blast']['blast']:
    # DATABASE: NCBI
	if config["blast"]["database"]=="ncbi":

		# Assign full NCBI taxonomic lineages to the BLAST hits using the
		# taxonomy IDs (staxid) and the lineage table (tax_lineage.h5).
		# Selected taxonomic classes can be dropped via the configuration.
		# Writes a temporary lineage table and a table containing all hits.
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
			log:
				os.path.join(config["general"]["output_dir"],"logfiles/blast/ncbi_taxonomy.log")
			conda:
				"../envs/classification/blast.yaml"
			script:
				"../scripts/classification/ncbi_taxonomy.py"

		# Combine annotated BLAST results with the selected abundance table
		# (SWARM, VSEARCH, or filtered input as selected below).
		# Write full, taxonomy-filtered, abundance, and metadata tables.
		# Taxonomy filtering is applied in the merge script.
		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv") \
                if config["clustering"] == "vsearch" and config["dataset"]["nanopore"] \
                else os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") \
                if config["general"]["seq_rep"] == "OTU" and not config["dataset"]["nanopore"] and config["clustering"] == "swarm" \
                else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_ncbi/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/merge_results.log")
			conda:
				"../envs/utilities/merge_results.yaml"
			script:
				"../scripts/utilities/merge_results_vsearch_blast.py" \
                if config["clustering"]== "vsearch" and config["dataset"]["nanopore"] \
                else "../scripts/utilities/merge_results.py"


    # DATABASE: SILVA
	elif config["blast"]["database"] == "silva":

		# Assign SILVA taxonomic lineages to the BLAST hits using the lineage
		# table (tax_lineage.h5). Selected taxonomic classes can be dropped
		# via the configuration. The result is a temporary lineage table.
		rule silva_taxonomy:
			input:
				blast_result = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv"),
				lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
			output:
				temp(os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv"))
			params:
				drop_tax_classes=str(config["blast"]["drop_tax_classes"])
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/silva_taxonomy.log")
			conda:
				"../envs/classification/blast.yaml"
			script:
				"../scripts/classification/silva_taxonomy.py"

		# Combine annotated BLAST results with the selected abundance table
		# (SWARM, VSEARCH, or filtered input as selected below).
		# Write full, taxonomy-filtered, abundance, and metadata tables.
		# Taxonomy filtering is applied in the merge script.
		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv") \
                if config["clustering"] == "vsearch" and config["dataset"]["nanopore"] \
                else os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") \
                if config["general"]["seq_rep"] == "OTU" and not config["dataset"]["nanopore"] and config["clustering"] == "swarm" \
                else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_silva/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/merge_results.log")
			conda:
				"../envs/utilities/merge_results.yaml"
			script:
				"../scripts/utilities/merge_results_vsearch_blast.py" \
                if config["clustering"]== "vsearch" and config["dataset"]["nanopore"] \
                else "../scripts/utilities/merge_results.py"


    # DATABASE: Eukaryome
	elif config["blast"]["database"] == "eukaryome":

		# Assign Eukaryome taxonomic lineages to the BLAST hits using the lineage
		# table (tax_lineage.h5). Selected taxonomic classes can be dropped
		# via the configuration. The result is a temporary lineage table.
		rule Eukaryome_taxonomy:
			input:
				blast_result = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv"),
				lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
			output:
				temp(os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv"))
			params:
				drop_tax_classes=str(config["blast"]["drop_tax_classes"])
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/Eukaryome_taxonomy.log")
			conda:
				"../envs/classification/blast.yaml"
			script:
				"../scripts/classification/eukaryome_taxonomy.py"

		# Combine annotated BLAST results with the selected abundance table
		# (SWARM, VSEARCH, or filtered input as selected below).
		# Write full, taxonomy-filtered, abundance, and metadata tables.
		# Taxonomy filtering is applied in the merge script.
		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv") \
                if config["clustering"] == "vsearch" and config["dataset"]["nanopore"] \
                else os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") \
                if config["general"]["seq_rep"] == "OTU" and not config["dataset"]["nanopore"] and config["clustering"] == "swarm" \
                else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_eukaryome/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_eukaryome/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_eukaryome/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_eukaryome/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/merge_results.log")
			conda:
				"../envs/utilities/merge_results.yaml"
			script:
				"../scripts/utilities/merge_results_vsearch_blast.py" \
                if config["clustering"]== "vsearch" and config["dataset"]["nanopore"] \
                else "../scripts/utilities/merge_results.py"


	# DATABASE: ROD
	elif config["blast"]["database"] == "rod":

		# Assign ROD taxonomic lineages to the BLAST hits using the lineage
		# table (tax_lineage.h5). The SILVA taxonomy script is reused here,
		# since both databases share the same lineage format.
		rule ROD_taxonomy:
			input:
				blast_result = os.path.join(config["general"]["output_dir"],"blast/blast_taxonomy.tsv"),
				lineage = os.path.join(os.path.dirname(config["blast"]["db_path"]), "tax_lineage.h5")
			output:
				temp(os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv"))
			params:
				drop_tax_classes=str(config["blast"]["drop_tax_classes"])
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/ROD_taxonomy.log")
			conda:
				"../envs/classification/blast.yaml"
			script:
				"../scripts/classification/silva_taxonomy.py"
		
        # Combine annotated BLAST results with the selected abundance table
		# (SWARM, VSEARCH, or filtered input as selected below).
		# Write full, taxonomy-filtered, abundance, and metadata tables.
		# Taxonomy filtering is applied in the merge script.
		rule merge_results:
			input:
				merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv") \
                if config["clustering"] == "vsearch" and config["dataset"]["nanopore"] \
                else os.path.join(config["general"]["output_dir"],"clustering/swarm_table.csv") \
                if config["general"]["seq_rep"] == "OTU" and not config["dataset"]["nanopore"] and config["clustering"] == "swarm" \
                else os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv"),
				blast_result=os.path.join(config["general"]["output_dir"],"blast/blast_taxonomic_lineage.tsv")
			output:
				complete=os.path.join(config["general"]["output_dir"],"finalData/blast_rod/full_table.csv"),
				filtered=os.path.join(config["general"]["output_dir"],"finalData/blast_rod/filtered_full_table.csv"),
				otus=os.path.join(config["general"]["output_dir"],"finalData/blast_rod/OTU_table.csv"),
				metadata=os.path.join(config["general"]["output_dir"],"finalData/blast_rod/metadata_table.csv"),
			params:
				seq_rep=str(config["general"]["seq_rep"]),
			log:
                os.path.join(config["general"]["output_dir"],"logfiles/blast/merge_results.log")
			conda:
				"../envs/utilities/merge_results.yaml"
			script:
				"../scripts/utilities/merge_results_vsearch_blast.py" \
                if config["clustering"]== "vsearch" and config["dataset"]["nanopore"] \
                else "../scripts/utilities/merge_results.py"
