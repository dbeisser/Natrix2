#!/bin/bash

# Check whether an argument has been passed
if [ -z "$1" ]; then
    varname="$PROJECT_NAME"
else
    varname="$1"
fi

# Check that varname is not empty
if [ -z "$varname" ]; then
    echo "Error: No project name provided."
    exit 1
fi

# Activate Conda environment natrix2
env_loc="$(conda info --base)/etc/profile.d/conda.sh"
source "$env_loc"
conda activate natrix2

# Check if the project name is 'test_docker' to modify the path
if [ "$varname" = "test_docker" ]; then
    yaml_file="docker_setup/test_docker.yaml"
else
    # Wait until the YAML file is found in the input directory
    while [ ! -f "input/$varname.yaml" ]; do
        echo "File input/$varname.yaml does not exist. Waiting 5 seconds."
        sleep 5
    done
    yaml_file="input/$varname.yaml"
fi

# Extract the number of cores from the YAML file
cores=$(python -c "import yaml; print(yaml.safe_load(open('$yaml_file'))['general']['cores'])")

# Abort if cores not found
if [ -z "$cores" ]; then
    echo "Error: No 'general.cores' entry found in $yaml_file"
    exit 1
fi

# Create the dataframe from samples
python create_dataframe.py "$yaml_file"

# Run Snakemake to start the analysis
snakemake --use-conda --cores $cores --configfile "$yaml_file"

# Keep container open
exec sh
