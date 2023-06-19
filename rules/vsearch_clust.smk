#Global # VSEARCH
rule vsearch_clust:
	input:
		fasta=os.path.join(config["general"]["output_dir"],"filtering/filtered.fasta"),
	output:
		os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus.fasta"),
		os.path.join(config["general"]["output_dir"],"clustering/vsearch_all_otus_tab.txt")
	threads: config["general"]["cores"]
	conda:
		"../envs/vsearch.yaml"
	shell:
		"""
			vsearch --cluster_size {input} \
			--threads {threads} \
			--id 0.97 \
			--strand plus \
			--sizein \
			--sizeout \
			--fasta_width 0 \
			--relabel OTU_ \
			--centroids {output[0]} \
			--otutabout {output[1]}
			--otutabout {output[1]}
		"""
