import pandas as pd
import collections as col

filtered_dict = pd.read_csv(snakemake.input["final_table_path2"], index_col="seqid").to_dict(orient="index")

with open(snakemake.input["merged"], "r") as swarms:
    seq_names = []
    for row in swarms:
        otu_seq_list = [s for s in row.split()]
        seq_names.append(otu_seq_list)

swarm_dict = col.defaultdict()
for j in range(len(seq_names)):
    seq_id = ">" + seq_names[j][0]
    if seq_id in filtered_dict:  # Check if sequence ID exists in filtered_dict
        swarm_dict[seq_id] = {}
        swarm_dict[seq_id]["sequences"] = filtered_dict[seq_id]["sequences"]
        for sample in filter(lambda i: i != "sequences", filtered_dict[seq_id].keys()):  # Iterate over samples
            swarm_dict[seq_id][sample] = sum([filtered_dict[">" + key][sample] for key in seq_names[j]])
        swarm_dict[seq_id]["seqid"] = "N{}_{}".format(seq_names[j][0].split(";")[0],
                                                      sum([value for key, value in swarm_dict[seq_id].items() if key not in {"seqid", "sequences"}]))

df = pd.DataFrame.from_dict(swarm_dict, orient="index").set_index("seqid")
df.to_csv(snakemake.output['out'])
