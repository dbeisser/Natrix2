import pathlib
import yaml
import pandas as pd
import numpy as np
from glob import glob
import os
import sys
import re

# Create the datatable containing the samples, units and paths of all
# fastq files formatted correctly. This is vital for the snakemake
# pipeline, without it, the wildcards can't be created.
# Additionally, options will be checked.

with open(sys.argv[1]) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)

def validate_sample_names(fl):
    # Regex for valid filenames: 'sampleX_(A|B)_R(1|2).fastq.gz'
    valid_pattern = re.compile(r'^[^_]+_(A|B)_R(1|2)\.fastq\.gz$')
    
    # Pre-check for disallowed patterns like "_AB_"
    invalid_samples = [
        file for file in fl
        if '_AB_' in file or not valid_pattern.match(file)
    ]
    
    if invalid_samples:
        error_message = (
            "\n\nERROR: Invalid Sample Names Detected!\n"
            + "\n".join(f"- {name}" for name in invalid_samples) +
            "\nREQUIREMENT: Sample names must follow the format:\n"
            "- Allowed: [sample_(A|B)_R(1|2).fastq.gz]\n"
            "- Not allowed: [_AB_], [_A_A_], [_A_B_], [_R3]\n"
        )
        raise ValueError(error_message)

def create_dataframe(fl, fpl, config, slice):
    validate_sample_names(fl) # Validate sample names before processing
    if config['merge']['paired_End'] and not config['general']['already_assembled']:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index =range(int(len(fl)/2)), dtype=str)
        i, j = (0, 0)

        while i < len(fl)/2:
            # last split needs to be fwd or rev read
            # second last can be unit
            unit = fl[j].split('_')[-2]
            if unit in ['A', 'B']:
                df.loc[i]['unit'] = unit
                df.loc[i]['sample'] = '_'.join(fl[j].split('_')[:-2])
            else:
                df.loc[i]['unit'] = ''
                df.loc[i]['sample'] = '_'.join(fl[j].split('_')[:-1])
            df.loc[i]['fq1'] = fpl[j][:slice]
            df.loc[i]['fq2'] = fpl[j+1][:slice]
            j += 2
            i += 1

    elif config['dataset']['nanopore']:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1'],
                          index=range(int(len(fl))), dtype=str)
        i = 0
        while i < len(fl):
            unit = fl[i].split('_')[-2]
            print(unit)
            # no units in sample name
            # print(fl[0])
            df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-2])
            df.loc[i]['fq1'] = fpl[i][:slice]
            df.loc[i]['unit'] = unit
            i += 1
    else:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index = range(int(len(fl))), dtype=str)

        i = 0
        while i < len(fl):
            # last split needs to be fwd or rev read
            # second last can be unit
            unit = fl[i].split('_')[-2]
            if unit in ['A', 'B']:
                df.loc[i]['unit'] = unit
                df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-2])
            else:
                df.loc[i]['unit'] = ''
                df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-1])
            df.loc[i]['fq1'] = fpl[i][:slice]
            df.loc[i]['fq2'] = np.nan
            i += 1
    return df


if __name__ == '__main__':
    # check config options
    if "-" in config["general"]["output_dir"]:
        sys.exit("Please rename output folder, do not use a dash in the folder name")

    # Check for uncompressed files
    uncompressed_files = glob(config['general']['filename'].rstrip("/") + '/*.fastq')
    if uncompressed_files:
        print("\nWarning: Uncompressed file(s) found!\n" +
            "\n".join(uncompressed_files) + 
            "\n\nPlease compress your file(s) first.\n" +
            "Command: {pigz -k filename.fastq}\n")
        sys.exit()

    if config["classify"]["mothur"] and config["blast"]["blast"]:
        sys.exit("Please decide whether to use blast or classification with mothur. Both config options cannot be set to TRUE")

    # Ensures [seq_rep: ASV] is not used with [postcluster: mumu]
    if config["general"]["seq_rep"] == "ASV" and config["postcluster"]["mumu"]:
        sys.exit("\n[Error]: Postclustering with [mumu] is not supported for [ASVs].\nPlease set [mumu] to FALSE in your config file. Workflow will be aborted.\n")

    # Ensures [seq_rep: ASV] is not used with [clustering: vsearch]
    if config['general']['seq_rep'] == 'ASV' and config['clustering'] == 'vsearch':
        sys.exit("\n[Error]: [seq_rep: ASV] cannot be used with [clustering: vsearch].\nPlease set [clustering] to [swarm] in your config file. Workflow will be aborted\n")

    if not config['general']['already_assembled']:
        file_path_list = [os.path.join(config["general"]["output_dir"],'demultiplexed/' + name.split('/')[-1]) for name in
                          sorted(glob(config['general']['filename'].rstrip("/") + '/*.gz'))]
        file_list = sorted([file_.split('/')[-1] for file_
                    in file_path_list])
        slice = -3 # Remove the .gz extension from the file paths.
    
    if config['dataset']['nanopore']:
        file_path_list = sorted(glob(os.path.join(config["general"]["filename"],'*R1.fastq.gz')))

        file_list = sorted([file_.split('/')[-1] for file_ in file_path_list])
        slice = None
        #print(file_list, file_path_list)

    # create dataframe
    df = create_dataframe(file_list, file_path_list, config, slice)
    print(df)

    pathlib.Path(config["general"]["output_dir"]).mkdir(parents=True, exist_ok=True)
    df.to_csv(os.path.join(config["general"]["output_dir"],config["general"]['units']), sep='\t')