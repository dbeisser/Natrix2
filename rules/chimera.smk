import os

# Detect potential chimeric sequences using VSEARCH (UCHIME3 de novo).
# Chimeras are artificial sequences that can arise during PCR amplification
# and may lead to an overestimation of microbial diversity.
# The rule performs reference-free chimera detection on the input sequences
# and separates them into chimeric and non-chimeric FASTA files.
# Additionally, it generates a UCHIME report with the detection results.
rule vsearch_chim:
    input:
        os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.dereplicated.fasta") if config["general"]["seq_rep"] == "OTU" else os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}_dada.fasta")
    output:
        uchime_out=temp(os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.vsearch.txt")),
        chim=os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}/{sample}_{unit}.chimera.fasta"),
        nonchim=os.path.join(config["general"]["output_dir"],"assembly/{sample}_{unit}.nonchimera.fasta")
    params:
        beta=config["chim"]["beta"],
        pseudo_c=config["chim"]["pseudo_count"],
        abskew=config["chim"]["abskew"]
    threads: 1
    log:
        os.path.join(config["general"]["output_dir"],"logfiles/chimera/vsearch_chim/{sample}_{unit}.log")
    conda:
         "../envs/analysis/vsearch.yaml"
    shell:
        """
            # Log stdout and stderr
            mkdir -p "$(dirname {log})"
            vsearch_log="{log}.vsearch"

            # Detect chimeras with VSEARCH
            vsearch --uchime3_denovo {input} \
                    -uchimeout {output.uchime_out} \
                    -chimeras {output.chim} \
                    -nonchimeras {output.nonchim} \
                    -xn {params.beta} \
                    -dn {params.pseudo_c} \
                    --threads {threads} \
                    -abskew {params.abskew} \
                    --log "$vsearch_log" 2>&1 | tee {log}
        """
