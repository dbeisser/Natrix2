#!/bin/bash

# Check whether an argument has been passed
if [ -z "$1" ]
then
    varname=$PROJECT_NAME
else
    varname=$1
fi

# Check that varname is not empty
if [ -z "$varname" ]
then
    echo "Error: No project name provided."
    exit 1
fi

env_loc=$(conda info --base)/etc/profile.d/conda.sh

# Activate Conda environment natrix
source $env_loc
conda activate natrix

# Wait until the YAML file is found in the input directory
while [ ! -f "input/$varname.yaml" ]
do
    echo "File input/$varname.yaml does not exist. Waiting 5 seconds."
    sleep 5
done

yaml_file="input/$varname.yaml"

# Extract the number of cores from the YAML file
cores=$(grep "cores:" $yaml_file | awk '{print $2}')

# Check whether the number of cores has been found
if [ -z "$cores" ]; then
    cores=1  # Fallback auf 1 Kern, wenn kein Wert gefunden wird
    echo "No 'cores' entry found in $yaml_file"
fi

# Create the dataframe with the Python script
python create_dataframe.py "$yaml_file"

# Run Snakemake to start the analysis
snakemake --use-conda --cores $cores --configfile "$yaml_file"

exec sh