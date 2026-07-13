#!/bin/bash

# Launcher
echo ""
echo "Natrix2 Pipeline Launcher"
echo "Exit launcher: enter 'exit' or 'quit'"
echo "Enter config:"
echo "Project:  your_config.yaml (root directory)"
echo "Test run: illumina_testrun (example config)"

# User input
read -e -p "> " varname

# Empty input
if [[ -z "$varname" ]]; then
    echo "ERROR: No configuration entered."
    exit 1
fi

# Exit launcher
case "${varname,,}" in
    exit|quit)
        echo "Launcher aborted."
        exit 0
        ;;
esac

# Append .yaml only if needed
if [[ "$varname" != *.yaml ]]; then
    varname="${varname}.yaml"
fi

# Built-in config presets
if [[ "$varname" == "illumina_testrun.yaml" ]]; then
    varname="config_presets/illumina_testrun.yaml"
fi

# Check if the config file exists
if [[ ! -f "$varname" ]]; then
    echo ""
    echo "ERROR: '$varname' not found. Please check the filename and try again."
    exit 1
fi

# Activate conda env
env_loc="$(conda info --base)/etc/profile.d/conda.sh"
source "$env_loc"
conda activate natrix2

# Extract number of cores from config file
cores=$(python -c "
import yaml
data = yaml.safe_load(open('$varname'))
print(data['general']['cores'])
" 2>/dev/null)

# Abort if cores not found
if [ -z "$cores" ]; then
    echo "ERROR: No 'general.cores' entry found in '$varname'."
    exit 1
fi

# Create dataframe and start pipeline in screen
if python create_dataframe.py "$varname"; then
    session_name="$(basename "$varname" .yaml)"
    screen -S "$session_name" bash -c "source \"$env_loc\"; conda activate natrix2; snakemake --use-conda --cores \"$cores\" --configfile \"$varname\" -p -r; exec sh"
else
    echo "ERROR: Failed to create dataframe."
    exit 1
fi
