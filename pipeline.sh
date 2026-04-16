#!/bin/bash

# Prints for pipeline.sh script
echo ""
echo "Natrix2 Pipeline Script"
echo "Enter project name (e.g. illumina_swarm):"

# Read user input file
read varname

# Append .yaml only if needed
if [[ "$varname" != *.yaml ]]; then
    varname="${varname}.yaml"
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
    session_name="${varname%.yaml}"
    screen -S "$session_name" bash -c "source \"$env_loc\"; conda activate natrix2; snakemake --use-conda --cores \"$cores\" --configfile \"$varname\" -p -r; exec sh"
fi
