import os

# Build and filter sequence-abundance tables from non-chimeric sample data.
# Optionally run AmpliconDuo analysis and prepare FASTA input for VSEARCH.

# Generate an unfiltered sequence-abundance table and an intermediate HDF5
# file from the non-chimeric FASTA files of the sequencing units.
rule unfiltered_table:
    input:
        expand(os.path.join(config["general"]["output_dir"],"assembly/{unit.sample}_{unit.unit}.nonchimera.fasta"), unit=units.reset_index().itertuples())
    output:
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_table.csv"),
        temp(os.path.join(config["general"]["output_dir"],"filtering/unfiltered_dict.hdf5"))
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/merging/unfiltered_table.log")
    conda:
        "../envs/utilities/unfiltered_table.yaml"
    script:
        "../scripts/utilities/unfiltered_table.py"

# Apply the configured filtering method and cutoff to the intermediate data.
# Write separate tables for retained and excluded sequences.
rule filtering:
    input:
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_dict.hdf5")
    output:
          temp(os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv")),
          os.path.join(config["general"]["output_dir"],"filtering/filtered_out_table.csv")
    params:
        filter_method=config["merge"]["filter_method"],
        cutoff=config["merge"]["cutoff"]
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/merging/filtering.log")
    conda:
        "../envs/classification/filtering.yaml"
    script:
        "../scripts/preprocessing/filtering.py"

# Run the AmpliconDuo R script using filtered and unfiltered tables.
# Pass plot settings and the multiple-testing correction method; save RData.
rule ampliconduo:
    input:
        filtered_table=os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv"),
        unfiltered_table=os.path.join(config["general"]["output_dir"],"filtering/unfiltered_table.csv")
    output:
        os.path.join(config["general"]["output_dir"],"filtering/figures/AmpliconDuo.RData")
    params:
        plot_ampduo=config["merge"]["plot_AmpDuo"],
        saving_format=config["merge"]["save_format"],
        p_corr=config["merge"]["ampli_corr"]
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/merging/ampliconduo.log")
    conda:
        "../envs/analysis/ampliconduo.yaml"
    script:
        "../scripts/analysis/ampliconduo.R"

# [VSEARCH-specific output preparation]
if config['clustering']=="vsearch":
    # Convert retained sequences to FASTA and add matching sequence IDs to
    # a CSV table. FASTA headers record each sequence's summed abundance.
    rule write_fasta:
        input:
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv")
        output:
            os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
            os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
        log:
            os.path.join(config["general"]["output_dir"], "logfiles/merging/write_fasta.log")
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
