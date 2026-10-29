import os

# Identify putative chimeric sequences with VSEARCH (UCHIME3 de novo).
# Chimeric sequences can arise during PCR amplification
# and may inflate estimated sequence diversity.
# The rule checks dereplicated OTU sequences or DADA-derived sequences
# and writes putative chimeric and non-chimeric sequences to FASTA files.
# It also writes a temporary UCHIME report.
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
            # Detect chimeras with VSEARCH; log stdout and stderr to file and terminal
            vsearch --uchime3_denovo {input} \
                --uchimeout {output.uchime_out} \
                --chimeras {output.chim} \
                --nonchimeras {output.nonchim} \
                --xn {params.beta} \
                --dn {params.pseudo_c} \
                --abskew {params.abskew} \
                --threads {threads} \
                2>&1 | tee {log:q}
        """
