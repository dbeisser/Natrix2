import os

rule unfiltered_table:
    input:
        expand(os.path.join(config["general"]["output_dir"],"assembly/{unit.sample}_{unit.unit}.nonchimera.fasta"), unit=units.reset_index().itertuples())
    output:
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_table.csv"),
        temp(os.path.join(config["general"]["output_dir"],"filtering/unfiltered_dict.hdf5"))
    conda:
        "../envs/unfiltered_table.yaml"
    script:
        "../scripts/unfiltered_table.py"

rule filtering:
    input:
        os.path.join(config["general"]["output_dir"],"filtering/unfiltered_dict.hdf5")
    output:
          temp(os.path.join(config["general"]["output_dir"],"filtering/filtered_table_temp.csv")),
          os.path.join(config["general"]["output_dir"],"filtering/filtered_out_table.csv")
    params:
        filter_method=config["merge"]["filter_method"],
        cutoff=config["merge"]["cutoff"]
    conda:
        "../envs/filtering.yaml"
    script:
        "../scripts/filtering.py"

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
    conda:
        "../envs/ampliconduo.yaml"
    log:
        os.path.join(config["general"]["output_dir"],"logs/ampliconduo.log")
    script:
        "../scripts/ampliconduo.R"

if config['clustering']=="vsearch":
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
