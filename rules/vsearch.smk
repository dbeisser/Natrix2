# Cluster filtered FASTA sequences with VSEARCH and write centroid,
# OTU-abundance, and cluster-assignment outputs.
rule vsearch_clust:
    input:
        fasta=os.path.join(config["general"]["output_dir"], "filtering/filtered.fasta"),
    output:
        os.path.join(config["general"]["output_dir"], "clustering/vsearch_all_otus.fasta"),
        os.path.join(config["general"]["output_dir"], "clustering/vsearch_all_otus_tab.txt"),
        os.path.join(config["general"]["output_dir"], "clustering/vsearch_uc.txt")
    threads: config["general"]["cores"]
    params:
        id=config["vsearch_id"]
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/vsearch/vsearch_clust.log")
    conda:
        "../envs/analysis/vsearch.yaml"
    shell:
        """
        # Cluster filtered sequences and write centroid, OTU abundance and assignment files.
        vsearch --cluster_fast {input.fasta:q} \
            --threads {threads} \
            --id {params.id} \
            --centroids {output[0]:q} \
            --otutabout {output[1]:q} \
            --uc {output[2]:q} 2>&1 | tee {log:q}
        """

# Prepare cluster identifiers for merging in the downstream script.
rule  write_clusters_uc:
    input:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus_tab.txt")
    output:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt")
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/vsearch/write_clusters_uc.log")
    script:
        "../scripts/analysis/vsearch_uc.py"

# Merge VSEARCH cluster identifiers with the filtered abundance table.
rule clust_merge_results_vsearch:
    input:
        merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt"),
        final_table_path2=os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
    output:
        out=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/vsearch/clust_merge_results_vsearch.log")
    conda:
        "../envs/utilities/merge_results.yaml"
    script:
        "../scripts/utilities/merge_clust_results.py"

# Generate a FASTA file with headers prepared for downstream analysis.
rule vsearch_header:
    input:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus.fasta"),
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
    output:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta")
    log:
        os.path.join(config["general"]["output_dir"], "logfiles/vsearch/vsearch_header.log")
    script:
        "../scripts/analysis/vsearch_fasta_headers.py"
