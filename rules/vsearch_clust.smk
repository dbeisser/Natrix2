#Global # VSEARCH
#rule combined_fasta:
#        input: 
rule vsearch_clust:
	input:
		fasta=os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
	output:
		os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus.fasta"),
		os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus_tab.txt"),
		os.path.join(config["general"]["output_dir"],"clustering/vsearch_uc.txt")
	threads: config["general"]["cores"]
	params:
	    id=config["vsearch_id"]
	conda:
		"../envs/vsearch.yaml"
	shell:
		"""
			vsearch --cluster_fast {input} \
			--threads {threads} \
			--id {params.id} \
			--centroids {output[0]} \
			--otutabout {output[1]}\
                        --uc {output[2]}
		"""

## to make cluster file for merging
rule  write_clusters_uc:
    input:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus_tab.txt")
    output:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt")
    script:
        "../scripts/vsearch_uc.py"

rule clust_merge_results_vsearch:
    input:
        merged=os.path.join(config["general"]["output_dir"],"clustering/vsearch_clusters_names.txt"),
        final_table_path2=os.path.join(config["general"]["output_dir"],"filtering/filtered_table.csv")
    output:
        out=os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
    conda:
        "../envs/merge_results.yaml"
    script:
        "../scripts/merge_clust_results.py"

rule vsearch_header:
    input:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus.fasta"),
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_table.csv")
    output:
        os.path.join(config["general"]["output_dir"],"clustering/vsearch_mod.fasta")
    script:
        "../scripts/vsearch_fasta_headers.py"
