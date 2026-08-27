import os 
sam=snakemake.input[0]
fasta=snakemake.input[1]
counts=snakemake.output[0]
rep_fasta=snakemake.output[1]

log1=os.system("samtools view -F 260 " + sam + " | cut -f 3 | sort | uniq -c | awk '{print $2,$1}' > " + counts)
log2=os.system("awk -F' ' '{ for(i=1; i<=$2;i++) print $1}' " + counts + " | awk '{print $1}' | xargs -n 1 -I '{}' grep -A 1 {} " + fasta + " > " + rep_fasta)

print(f"counts generated\n")
print(f"repeat fasta generated\n")
