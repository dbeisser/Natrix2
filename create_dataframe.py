import pathlib
import yaml
import pandas as pd
import numpy as np
from glob import glob
import os
import sys
import re

# Build the sample table containing sample IDs, sequencing units, and FASTQ paths.
# Snakemake uses this table to resolve wildcards and locate the input files.
# The script also validates relevant configuration options and input metadata.

with open(sys.argv[1]) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)

def validate_sample_names(fl):
    # Accept only the documented format: sampleID_(A|B)_R(1|2).fastq.gz.
    valid_pattern = re.compile(r'^[^_]+_(A|B)_R(1|2)\.fastq\.gz$')
    
    # Collect filenames that either contain a forbidden unit pattern or do not
    # match the required naming convention.
    invalid_samples = [
        file for file in fl
        if '_AB_' in file or not valid_pattern.match(file)
    ]
    
    if invalid_samples:
        error_message = (
            "\nERROR: Invalid FASTQ sample names detected.\n"
            + "\n".join(f"- {name}" for name in invalid_samples) +
            "\nRequired filename format:\n"
            "- Allowed: sample_(A|B)_R(1|2).fastq.gz\n"
            "- Not allowed: _AB_, _A_A_, _A_B_, or _R3\n"
        )
        raise ValueError(error_message)


def validate_primer_samples(file_list, config):
    # Read the primer-table path from the general configuration section.
    try:
        primer_table_path = config["general"]["primertable"]
    except KeyError:
        sys.exit(
            "\nERROR: Primer-table path is missing in the config.\n"
            "Expected config entry:\n"
            "general:\n"
            "  primertable: path/to/primer_table.csv\n"
        )

    # Reject missing, empty, or non-string primer-table paths before file access.
    if not isinstance(primer_table_path, str) or not primer_table_path.strip():
        sys.exit(
            "\nERROR: Config entry 'general: primertable' is empty or invalid.\n"
        )

    # Confirm that the configured primer table exists.
    if not os.path.isfile(primer_table_path):
        sys.exit(
            "\nERROR: Primer table was not found:\n"
            f"- {primer_table_path}\n"
        )

    try:
        # Let pandas detect comma-, semicolon-, or tab-separated input.
        primer_df = pd.read_csv(
            primer_table_path,
            sep=None,
            engine="python",
            dtype=str
        )
    except Exception as exc:
        sys.exit(
            "\nERROR: Primer table could not be read:\n"
            f"- {primer_table_path}\n"
            f"Reason: {exc}\n"
        )

    # Normalize column names by removing surrounding whitespace and a possible
    # UTF-8 byte-order mark.
    primer_df.columns = [
        str(column).strip().lstrip("\ufeff")
        for column in primer_df.columns
    ]

    if "Probe" not in primer_df.columns:
        sys.exit(
            "\nERROR: Primer table must contain a column named 'Probe'.\n"
            f"Primer table: {primer_table_path}\n"
            f"Available columns: {', '.join(primer_df.columns)}\n"
        )

    # Map each FASTQ filename to its normalized sample name without the read suffix.
    file_to_sample = {}
    invalid_input_files = []

    for filename in file_list:
        basename = os.path.basename(filename)
        match = re.fullmatch(r"(?P<sample>.+)_R[12]\.fastq\.gz", basename)

        if match is None:
            invalid_input_files.append(basename)
        else:
            file_to_sample[basename] = match.group("sample")

    if invalid_input_files:
        sys.exit(
            "\nERROR: Sample names could not be extracted from these files:\n- "
            + "\n- ".join(sorted(invalid_input_files))
            + "\nExpected filename ending: _R1.fastq.gz or _R2.fastq.gz\n"
        )

    # Keep both list and set representations for duplicate and membership checks.
    input_sample_list = list(file_to_sample.values())
    input_samples = set(input_sample_list)

    # Normalize Probe values before comparing them with FASTQ-derived names.
    primer_series = (
        primer_df["Probe"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    empty_rows = primer_series[primer_series == ""].index.tolist()
    if empty_rows:
        sys.exit(
            "\nERROR: Empty entries were found in primer-table column 'Probe'.\n"
            "Affected CSV row(s): "
            + ", ".join(str(row_index + 2) for row_index in empty_rows)
            + "\n"
        )

    primer_sample_list = primer_series.tolist()
    primer_samples = set(primer_sample_list)

    # Primer-table entries must describe samples, not individual read files.
    invalid_primer_names = sorted(
        name for name in primer_samples
        if re.search(r"_R[12](?:\.fastq(?:\.gz)?)?$", name)
    )
    if invalid_primer_names:
        sys.exit(
            "\nERROR: Primer-table names must not contain _R1 or _R2.\n"
            "Correct example: Probe1_A\n"
            "Incorrect example: Probe1_A_R1\n"
            "Invalid primer-table entries:\n- "
            + "\n- ".join(invalid_primer_names)
            + "\n"
        )

    duplicate_input_samples = sorted({
        sample for sample in input_sample_list
        if input_sample_list.count(sample) > 1
    })

    # In paired-end datasets, R1 and R2 intentionally map to the same normalized
    # sample name. Duplicate normalized names are therefore relevant only for
    # datasets that are not configured as paired-end.
    if config["merge"]["paired_End"]:
        duplicate_input_samples = []

    # Detect repeated Probe entries because each sample should occur only once.
    duplicate_primer_samples = sorted({
        sample for sample in primer_sample_list
        if primer_sample_list.count(sample) > 1
    })

    # Compare both sources in both directions to report all mismatches.
    missing_in_primer = sorted(input_samples - primer_samples)
    missing_in_input = sorted(primer_samples - input_samples)

    # Collect all consistency problems so the user receives one complete report.
    errors = []

    if duplicate_input_samples:
        errors.append(
            "Duplicate normalized sample names in input files:\n- "
            + "\n- ".join(duplicate_input_samples)
        )

    if duplicate_primer_samples:
        errors.append(
            "Duplicate sample names in primer table:\n- "
            + "\n- ".join(duplicate_primer_samples)
        )

    if missing_in_primer:
        affected_files = sorted(
            filename
            for filename, sample in file_to_sample.items()
            if sample in missing_in_primer
        )
        errors.append(
            "Input samples missing from primer table:\n- "
            + "\n- ".join(missing_in_primer)
            + "\nAffected input files:\n- "
            + "\n- ".join(affected_files)
        )

    if missing_in_input:
        errors.append(
            "Primer-table samples without matching input files:\n- "
            + "\n- ".join(missing_in_input)
        )

    # Stop before writing units.tsv when any sample inconsistency is found.
    if errors:
        sys.exit(
            "\nERROR: Sample-name validation failed.\n\n"
            + "\n\n".join(errors)
            + f"\n\n-Unique input samples: {len(input_samples)}"
            + f"\n-Primer-table samples: {len(primer_samples)}"
            + "\n\nunits.tsv file was NOT created.\n"
        )

    # Report the validated sample counts before units.tsv is created.
    print(
        "\nSample-name validation successful.\n"
        f"-Input FASTQ files: {len(file_list)}\n"
        f"-Unique input samples: {len(input_samples)}\n"
        f"-Primer-table samples: {len(primer_samples)}\n"
        "All sample names match exactly.\n"
        "units.tsv file will now be created.\n"
    )

def create_dataframe(fl, fpl, config, slice):
    # Validate all FASTQ filenames before constructing the sample table.
    validate_sample_names(fl)
    if config['merge']['paired_End'] and not config['general']['already_assembled']:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index =range(int(len(fl)/2)), dtype=str)
        i, j = (0, 0)

        while i < len(fl)/2:
            # The final filename component identifies the read direction, while
            # the preceding component may identify sequencing unit A or B.
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
            # Nanopore filenames are parsed using the same documented unit field.
            df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-2])
            df.loc[i]['fq1'] = fpl[i][:slice]
            df.loc[i]['unit'] = unit
            i += 1
    else:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index = range(int(len(fl))), dtype=str)

        i = 0
        while i < len(fl):
            # The final filename component identifies the read direction, while
            # the preceding component may identify sequencing unit A or B.
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
    # Validate configuration combinations before processing input files.
    if "-" in config["general"]["output_dir"]:
        sys.exit(
            "\nERROR: The output directory name contains a dash.\n"
            "Rename the output directory and use a name without dashes.\n"
        )

    # Stop when uncompressed FASTQ files are present because the workflow
    # expects gzip-compressed input.
    uncompressed_files = glob(config['general']['filename'].rstrip("/") + '/*.fastq')
    if uncompressed_files:
        print(
            "\nERROR: Uncompressed FASTQ files were found.\n"
            + "\n".join(f"- {file}" for file in uncompressed_files)
            + "\n\nCompress these files before starting the workflow.\n"
            + "Command: pigz -k filename.fastq\n"
        )
        sys.exit()

    if config["classify"]["mothur"] and config["blast"]["blast"]:
        sys.exit(
            "\nERROR: BLAST and Mothur classification cannot both be enabled.\n"
            "Set either 'blast: blast' or 'classify: mothur' to FALSE.\n"
        )

    # ASV output is incompatible with MUMU post-clustering.
    if config["general"]["seq_rep"] == "ASV" and config["postcluster"]["mumu"]:
        sys.exit(
            "\nERROR: MUMU post-clustering is not supported for ASVs.\n"
            "Set 'postcluster: mumu' to FALSE in the configuration file.\n"
        )

    # ASV output is incompatible with VSEARCH clustering in this workflow.
    if config['general']['seq_rep'] == 'ASV' and config['clustering'] == 'vsearch':
        sys.exit(
            "\nERROR: ASV representation cannot be used with VSEARCH clustering.\n"
            "Set 'clustering' to 'swarm' in the configuration file.\n"
        )

    if not config['general']['already_assembled']:
        file_path_list = [os.path.join(config["general"]["output_dir"],'demultiplexed/' + name.split('/')[-1]) for name in
                          sorted(glob(config['general']['filename'].rstrip("/") + '/*.gz'))]
        file_list = sorted([file_.split('/')[-1] for file_
                    in file_path_list])
        # Remove the .gz suffix from paths passed to downstream rules.
        slice = -3
    
    if config['dataset']['nanopore']:
        file_path_list = sorted(glob(os.path.join(config["general"]["filename"],'*R1.fastq.gz')))

        file_list = sorted([file_.split('/')[-1] for file_ in file_path_list])
        slice = None
        # Keep complete paths for Nanopore input files.

    # Create the sample table used by Snakemake.
    df = create_dataframe(file_list, file_path_list, config, slice)
    print(df)

    # Verify that normalized FASTQ sample names and primer-table Probe entries
    # match exactly before writing units.tsv.
    validate_primer_samples(file_list, config)

    pathlib.Path(config["general"]["output_dir"]).mkdir(parents=True, exist_ok=True)
    df.to_csv(os.path.join(config["general"]["output_dir"],config["general"]['units']), sep='\t')
