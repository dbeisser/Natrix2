import csv

def update_fasta_headers(fasta_file, csv_file):
    # Read the CSV file and create a dictionary mapping the new headers
    header_mapping = {}
    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            mod=row[0].split('_')
            modified_header = mod[0].replace("N",">")
            header_mapping[modified_header] = row[0]
    # Open the FASTA file and create new fasta file with header same as vsearch_table.csv
    temp_file = open(snakemake.output[0], "w") #new file
    with open(fasta_file_path, 'r') as fasta, temp_file as temp:
        for line in fasta:
            if line.startswith('>'):
                header=line.strip()
                header1=header.split(';')
 #               print(header1[0])
                if header1[0] in header_mapping:
                    modified_header = header_mapping[header1[0]]
                    temp.write(f'>{modified_header}\n')
                else:
                    temp.write(line)
            else:
                temp.write(line)


    print('FASTA headers have been updated successfully.')

# Provide the paths to your FASTA file and CSV file
fasta_file_path = snakemake.input[0]
csv_file_path = snakemake.input[1]

# Call the function to update the headers
update_fasta_headers(fasta_file_path, csv_file_path)
