def parse_uc_file(uc_file):
    clusters = {}
    with open(uc_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if fields[0] == 'S':
                cluster_id = fields[8]
            if cluster_id != fields[8]:
                sequence_id = fields[8]
            else:
                sequence_id=''

            clusters.setdefault(cluster_id, []).append(sequence_id)
    return clusters

def write_cluster_file(clusters, output_file):
    with open(output_file, 'w') as file:
        for cluster_id, sequences in clusters.items():
            file.write(cluster_id + '\t' + ' '.join(sequences) + '\n')

uc_file = snakemake.input[0] # Replace with the actual UC file path
output_file = snakemake.output[0] # Replace with the desired output file path

clusters = parse_uc_file(uc_file)
write_cluster_file(clusters, output_file)
