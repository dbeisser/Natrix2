#!/bin/bash

# Prints for pipeline.sh script
echo ""
echo "Natrix2 Pipeline Script"
echo "Enter project name (e.g. illumina_swarm):"

# Read user input file
read varname

# Activate conda env
env_loc=$(conda info --base)/etc/profile.d/conda.sh
source $env_loc
conda activate natrix

# Append .yaml to the project name
varname+=".yaml"

# Check if the config file exists
if [[ ! -f "$varname" ]]; then
    echo ""
    echo "ERROR: '$varname' not found. Please check the filename and try again."
    exit 1
fi

# Extract number of cores from config file
cores=$(grep "cores : " $varname | awk '{print $3}')

# Run Python pre-step and start pipeline in screen
if python create_dataframe.py "$varname"; then
    screen -S $varname bash -c "source $env_loc;conda activate natrix;snakemake --use-conda --cores $cores --configfile $varname -p -r  ; exec sh"
fi
