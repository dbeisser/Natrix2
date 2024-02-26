import pandas as pd
import re
import sys

otu_file=open(snakemake.input[1], "r")
taxonomy=open(snakemake.input[0], "r")

output_all=open(snakemake.output[0], "w")
output_table=open(snakemake.output[1], "w")
output_meta=open(snakemake.output[2], "w")


taxonomy = pd.read_csv(taxonomy, sep='\t', names=['seqid', 'taxonomy']) # read mothur taxonomy file
if snakemake.params['clustering'] == 'swarm' or snakemake.params['clustering'] == 'vsearch':
    taxonomy['seqid'] = '>' + taxonomy['seqid'].astype(str)
#print(taxonomy.head())
otu_count = pd.read_csv(otu_file) #read otu table

tax=taxonomy['taxonomy']
ID=taxonomy['seqid']
concat_col=pd.concat([ID, tax], axis=1) #concatenate the id and taxonomy column

otu_count['seqid']=otu_count['seqid'] # change ID if ASV table

#for item in otu_count['seqid']:
#	print(item)


final_file=pd.merge(otu_count, concat_col, left_on='seqid', right_on='seqid', how='left') # merge mothur taxonomy to otu abundance file

final_file.to_csv(output_all, index=False) #output file

out2 = final_file[["seqid", "sequences", "taxonomy"]]

# substitute bootstrap values if present and append to dataframe
taxonomy2 = out2["taxonomy"].tolist()
if  isinstance(taxonomy2[0], str) and re.match("(\d*)", taxonomy2[0]):
    taxonomy2 = list(map(lambda x: x if pd.isna(x) else re.sub(r"\([^()]*\)", "", x), taxonomy2))
    pd.options.mode.chained_assignment = None
    out2.loc[:,'taxonomy2'] = taxonomy2

out2.to_csv(output_meta, index=False)

# get only counts
out3=otu_count.drop(["sequences"], axis=1)
out3.to_csv(output_table, index=False)
