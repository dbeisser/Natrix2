import os
import subprocess
import glob
import sys

def check_seqkit_installed():
    try:
        result = subprocess.run(["seqkit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print("ERROR: seqkit is not installed. Please install seqkit first.")
            sys.exit(1)
    except FileNotFoundError:
        print("ERROR: seqkit is not installed. Please install seqkit first.")
        sys.exit(1)

def get_seqkit_sequence_count(file_path):
    result = subprocess.run(f"seqkit stats {file_path} | awk 'NR==2 {{print $4}}'", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True)
    if result.returncode != 0:
        print(f"ERROR: Processing file {file_path}: {result.stdout}")
        return 0
    try:
        return int(result.stdout.strip().replace(',', ''))
    except ValueError as e:
        print(f"ERROR: Parsing output for {file_path}: {result.stdout} - {e}")
        return 0

def check_sequences_in_folder(folder_path, threshold=150):
    fastq_files = glob.glob(os.path.join(folder_path, '*.fastq')) + glob.glob(os.path.join(folder_path, '*.fastq.gz'))

    if not fastq_files:
        print(f"\nERROR:\nNo files (fastq or fastq.gz) found in the directory.\nfolder_path: {folder_path}\n")
        return

    files_below_threshold = []

    print("\nANALYZE: files(*.fastq, *.fastq.gz)...")
    for file_path in fastq_files:
        sequence_count = get_seqkit_sequence_count(file_path)
        print(f"{os.path.basename(file_path)}: {sequence_count} sequences")
        if sequence_count < threshold:
            files_below_threshold.append((os.path.basename(file_path), sequence_count))

    if files_below_threshold:
        print("\nWARNING: The following files have sequence counts below the threshold of {}:".format(threshold))
        for file_name, count in files_below_threshold:
            print(f"{file_name}: {count} sequences. Please review and consider removing it.")
    else:
        print(f"\nNOTE: All files have sequence counts above the threshold of {threshold}.")
    print("")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        script_name = os.path.basename(__file__)
        print(f"Usage: python3 {script_name} <folder_path> <threshold>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    threshold = int(sys.argv[2])

    check_seqkit_installed()
    check_sequences_in_folder(folder_path, threshold)
