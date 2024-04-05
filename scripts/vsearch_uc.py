def write_clusters_uc(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    sample_ids = lines[0].strip().split('\t')[1:]

    with open(output_file, 'w') as f:
        # Iterate through lines starting from the second line
        for line in lines[1:]:
            elements = line.strip().split('\t')
            otu_id = elements[0]
            abundances = elements[1:]
            # Create lists to store sample IDs and abundances
            samples = []
            non_matching_samples = []
            # Iterate through abundances
            for idx, abundance in enumerate(abundances):
                # If abundance is greater than 0
                if int(abundance) > 0:
                    sample_id = sample_ids[idx]
                    if sample_id == otu_id:
                        samples.append((sample_id, abundance))
                    else:
                        non_matching_samples.append((sample_id, abundance))
            # Sort the samples based on the cluster ID
            sorted_samples = sorted(samples, key=lambda x: x[0])
            sorted_non_matching_samples = sorted(non_matching_samples, key=lambda x: x[0])
            # Construct the output line
            output_line = '; '.join([f'{sample[0]};size={sample[1]}' for sample in sorted_samples])
            non_matching_output_line = '; '.join([f'{sample[0]};size={sample[1]}' for sample in sorted_non_matching_samples])
            # Write the output line
            if non_matching_output_line:
                if output_line:
                    f.write(f'{output_line}; {non_matching_output_line};\n')
                else:
                    f.write(f'{non_matching_output_line};\n')
            elif output_line:
                f.write(f'{output_line};\n')

input_file = snakemake.input[0]
output_file = snakemake.output[0]
write_clusters_uc(input_file, output_file)
