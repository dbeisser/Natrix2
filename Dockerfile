# Sets base image to Miniconda3
FROM continuumio/miniconda3
# Changes default shell to Bash
SHELL ["/bin/bash", "-c"]
# Copies current directory to /app
COPY . /app
# Sets working directory to /app
WORKDIR /app
# Update and clean the system packages
RUN apt update && apt-get install -y libltdl7 && apt upgrade -y && apt-get purge -y && apt-get clean
# Create the environment
RUN conda env create -f natrix.yaml
# Makes the docker_pipeline.sh script executable
RUN chmod +x docker_pipeline.sh

# Execute docker_dummyfiles processing workflows and remove temporary files
RUN env_loc=$(conda info --base)/etc/profile.d/conda.sh && \
    source $env_loc && \ 
    conda activate natrix && \
    mkdir docker_dummy_env1 && \
    touch docker_dummy_env1.csv && \
    cp docker_dummyfiles/units.tsv docker_dummy.tsv && \
    mkdir docker_dummy_env2 && \
    touch docker_dummy_env2.csv && \
    python create_dataframe.py docker_dummyfiles/docker_dummy_env1.yaml && \
    snakemake --configfile docker_dummyfiles/docker_dummy_env1.yaml --cores 1 --use-conda --conda-create-envs-only && \
    python create_dataframe.py docker_dummyfiles/docker_dummy_env2.yaml && \
    snakemake --configfile docker_dummyfiles/docker_dummy_env2.yaml --cores 1 --use-conda --conda-create-envs-only && \
    rm -rf docker_dummy_env1 && \
    rm docker_dummy_env1.csv && \
    rm -rf docker_dummy_env2 && \
    rm docker_dummy_env2.csv && \
    rm docker_dummy.tsv && \
    rm -rf Illumina_results_swarm

# Run the pipeline script with the specified project name
CMD ["sh","-c", "./docker_pipeline.sh $PROJECT_NAME"]

# Ensures all RUN commands use the 'myenv' Conda environment. Currently not needed.
# SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]