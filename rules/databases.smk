import os

# Download and prepare reference databases for the configured classifier.
# Only the branch matching classify.database is defined by this file.
if config["classify"]["database"] == "pr2":
    # Download and decompress PR2 reference sequences and taxonomy files.
    rule download_pr2:
        output:
            expand(["database/pr2db.{pr2_db_version}.fasta", "database/pr2db.{pr2_db_version}.tax"], pr2_db_version=config["database_version"]["pr2"])
        params:
            db_version=config["database_version"]["pr2"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/databases/download_pr2.log")
        shell:
            """
                # Logging: stdout and stderr of each command are written to the log file and terminal;
                # the first command creates the log (tee), all following commands append (tee -a).
                # For pigz/gunzip/cat the data still go to the output file, only messages are logged.
                # Download compressed PR2 reference sequences into the working directory.
                wget -P ./ --progress=bar https://github.com/pr2database/pr2database/releases/download/v{params.db_version}/pr2_version_{params.db_version}_SSU_mothur.fasta.gz 2>&1 | tee {log:q};
                # Download the matching compressed PR2 taxonomy file.
                wget -P ./ --progress=bar https://github.com/pr2database/pr2database/releases/download/v{params.db_version}/pr2_version_{params.db_version}_SSU_mothur.tax.gz 2>&1 | tee -a {log:q};
                # Decompress the FASTA into the database directory.
                pigz -dc ./pr2_version_{params.db_version}_SSU_mothur.fasta.gz 2>&1 > database/pr2db.{params.db_version}.fasta | tee -a {log:q};
                # Decompress the taxonomy file into the database directory.
                pigz -dc ./pr2_version_{params.db_version}_SSU_mothur.tax.gz 2>&1 > database/pr2db.{params.db_version}.tax | tee -a {log:q};
                # Remove the downloaded compressed archives.
                rm ./pr2_version_{params.db_version}_SSU_mothur.tax.gz ./pr2_version_{params.db_version}_SSU_mothur.fasta.gz 2>&1 | tee -a {log:q};
            """

elif config["classify"]["database"] == "unite":
    # Download and extract the specified UNITE release in mothur format.
    rule download_unite:
        output:
            "database/unite_v10.fasta",
            "database/unite_v10.tax"
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/databases/download_unite.log")
        shell:
            """
                # Logging: stdout and stderr of each command are written to the log file and terminal;
                # the first command creates the log (tee), all following commands append (tee -a).
                # For pigz/gunzip/cat the data still go to the output file, only messages are logged.
                # Download the UNITE mothur release archive.
                wget -P ./ --progress=bar  -O sh_mothur_release_2024-04-04.tgz --progress=bar https://s3.hpc.ut.ee/plutof-public/original/10220f9f-17c8-4a2f-aada-1f852f50c0f7.tgz 2>&1 | tee {log:q};
                # Create the directory used to unpack the release.
                mkdir -p sh_mothur_release_2024-04-04 2>&1 | tee -a {log:q};
                # Extract the archive into that directory.
                tar -xvf sh_mothur_release_2024-04-04.tgz -C sh_mothur_release_2024-04-04 2>&1 | tee -a {log:q};
                # Move the taxonomy mapping to its expected output path.
                mv sh_mothur_release_2024-04-04/UNITEv10_sh_99.tax database/unite_v10.tax 2>&1 | tee -a {log:q};
                # Move the reference FASTA to its expected output path.
                mv sh_mothur_release_2024-04-04/UNITEv10_sh_99.fasta {output[0]} 2>&1 | tee -a {log:q};
                # Remove the extracted release directory.
                rm -rf sh_mothur_release_2024-04-04 2>&1 | tee -a {log:q};
                # Remove the downloaded archive.
                rm sh_mothur_release_2024-04-04.tgz 2>&1 | tee -a {log:q};
            """
            # IMPORTANT: By using this direct download link (https://s3.hpc.ut.ee/plutof-public/original/10220f9f-17c8-4a2f-aada-1f852f50c0f7.tgz), you agree to the terms and conditions of the UNITE website.
            # Please review the terms and conditions at: https://doi.plutof.ut.ee/doi/10.15156/BIO/2959342
            # License: Data is available under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
            # Reference: Abarenkov K, Nilsson RH, Larsson K-H, Taylor AFS, May TW, Frøslev TG, Pawlowska J, Lindahl B, Põldmaa K, Truong C, Vu D, Hosoya T, Niskanen T, Piirmann T, Ivanov F, Zirk A, Peterson M, 
            # Cheeke TE, Ishigami Y, Jansson AT, Jeppesen TS, Kristiansson E, Mikryukov V, Miller JT, Oono R, Ossandon FJ, Paupério J, Saar I, Schigel D, Suija A, Tedersoo L, Kõljalg U. 2023. 
            # The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Research, https://doi.org/10.1093/nar/gkad1039

elif config["classify"]["database"] == "silva":
        # Download and decompress SILVA sequences and taxonomy mapping data.
        rule download_silva:
            output:
                expand(["database/silva_db.{silva_db_version}.fasta", "database/silva_db.{silva_db_version}.tax.temp"], silva_db_version=config["database_version"]["silva"])
            params:
                db_version=config["database_version"]["silva"]
            log:
                os.path.join(config["general"]["output_dir"], "logfiles/databases/download_silva.log")
            shell:
                """
                # Logging: stdout and stderr of each command are written to the log file and terminal;
                # the first command creates the log (tee), all following commands append (tee -a).
                # For pigz/gunzip/cat the data still go to the output file, only messages are logged.
                # Download the compressed SILVA reference FASTA.
                wget -P ./ --progress=bar -O database/silva_{params.db_version}.fasta.gz https://www.arb-silva.de/fileadmin/silva_databases/release_{params.db_version}/Exports/SILVA_{params.db_version}_SSURef_tax_silva.fasta.gz 2>&1 | tee {log:q};
                # Decompress the reference sequences into the expected FASTA.
                pigz -dc database/silva_{params.db_version}.fasta.gz 2>&1 > database/silva_db.{params.db_version}.fasta | tee -a {log:q};
                # Download the SILVA taxonomy mapping archive.
                wget -P ./ --progress=bar -O database/silva_{params.db_version}.tax.gz https://www.arb-silva.de/fileadmin/silva_databases/release_{params.db_version}/Exports/taxonomy/taxmap_slv_ssu_ref_{params.db_version}.txt.gz 2>&1 | tee -a {log:q}
                # Decompress the taxonomy mapping into an intermediate file.
                pigz -dc database/silva_{params.db_version}.tax.gz 2>&1 > database/silva_db.{params.db_version}.tax.temp | tee -a {log:q};
                """
        
        # Convert the downloaded SILVA taxonomy mapping using a Python script.
        rule edit_silva:
            input:
                expand("database/silva_db.{silva_db_version}.tax.temp", silva_db_version=config["database_version"]["silva"])
            output:
                expand("database/silva_db.{silva_db_version}.tax", silva_db_version=config["database_version"]["silva"])
            log:
                os.path.join(config["general"]["output_dir"], "logfiles/databases/edit_silva.log")
            conda:
                "../envs/classification/blast.yaml"
            script:
                "../scripts/classification/edit_silva_mothur.py"

elif config["classify"]["database"] == "rod":
    # Download and extract ROD reference sequences and tabular annotations.
    rule download_ROD:
        output:
            expand(["database/RODdb.{ROD_db_version}_reference_sequences.fasta"], ROD_db_version=config["database_version"]["rod"]), 
            temp(expand(["database/RODdb.{ROD_db_version}_reference_sequences.tab.temp"], ROD_db_version=config["database_version"]["rod"]))
        params:
            db_version=config["database_version"]["rod"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/databases/download_ROD.log")
        shell:
            """
                # Logging: stdout and stderr of each command are written to the log file and terminal;
                # the first command creates the log (tee), all following commands append (tee -a).
                # For pigz/gunzip/cat the data still go to the output file, only messages are logged.
                # Download the compressed ROD release archive.
                wget -P ./ --progress=bar https://github.com/krabberod/ROD/archive/refs/tags/v{params.db_version}-goldenrod.zip 2>&1 | tee {log:q};
                # Extract the release into the database directory.
                unzip -o ./v{params.db_version}-goldenrod.zip -d ./database/ 2>&1 | tee -a {log:q};
                # Decompress the ROD reference sequences into the expected FASTA.
                gunzip -c database/ROD-{params.db_version}-goldenrod/ROD_v{params.db_version}_reference_sequences.fasta.gz 2>&1 > database/RODdb.{params.db_version}_reference_sequences.fasta | tee -a {log:q}
                # Decompress the ROD annotation table into a temporary file.
                gunzip -c database/ROD-{params.db_version}-goldenrod/ROD_v{params.db_version}_reference_sequences.tab.gz 2>&1 > database/RODdb.{params.db_version}_reference_sequences.tab.temp | tee -a {log:q}
                # Remove the downloaded ZIP archive.
                rm ./v{params.db_version}-goldenrod.zip 2>&1 | tee -a {log:q};
                # Remove the extracted release directory.
                rm -rf database/ROD-{params.db_version}-goldenrod 2>&1 | tee -a {log:q}
            """

    # Convert the intermediate ROD annotation table using a Python script.
    rule edit_ROD:
        input:
            expand("database/RODdb.{ROD_db_version}_reference_sequences.tab.temp", ROD_db_version=config["database_version"]["rod"])
        output:
            expand("database/RODdb.{ROD_db_version}_reference_sequences.tax", ROD_db_version=config["database_version"]["rod"])
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/databases/edit_ROD.log")
        conda:
            "../envs/classification/blast.yaml"
        script:
            "../scripts/classification/edit_ROD_mothur.py"

elif config["classify"]["database"] == "eukaryome":
    # Download and extract Eukaryome mothur reference and taxonomy files.
    rule download_eukaryome:
        output:
            expand(["database/EUK.{eukaryome_db_version}.tax", "database/EUK.{eukaryome_db_version}.fasta"], eukaryome_db_version=config["database_version"]["eukaryome"])
        params:
            db_version=config["database_version"]["eukaryome"],
            eukaryome_version=config["database_path"]["eukaryome_version"]
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/databases/download_eukaryome.log")
        shell:
            """
                # Logging: stdout and stderr of each command are written to the log file and terminal;
                # the first command creates the log (tee), all following commands append (tee -a).
                # For pigz/gunzip/cat the data still go to the output file, only messages are logged.
                # Download the selected Eukaryome mothur archive.
                wget -P ./database/ --progress=bar https://sisu.ut.ee/wp-content/uploads/sites/643/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip 2>&1 | tee {log:q}
                # Extract the ZIP archive into the database directory.
                python3 -m zipfile -e database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip database/ 2>&1 | tee -a {log:q}
                # Copy the extracted FASTA to the configured output filename.
                cat database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.fasta 2>&1 > database/EUK.{params.db_version}.fasta | tee -a {log:q}
                # Copy the extracted taxonomy file to the configured output filename.
                cat database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.tax 2>&1 > database/EUK.{params.db_version}.tax | tee -a {log:q}
                # Delete the downloaded ZIP archive.
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.zip 2>&1 | tee -a {log:q}
                # Delete the extracted FASTA after copying.
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.fasta 2>&1 | tee -a {log:q}
                # Delete the extracted taxonomy file after copying.
                rm database/mothur_EUK_{params.eukaryome_version}_v{params.db_version}.tax 2>&1 | tee -a {log:q}
            """
