# Base image with Miniconda
FROM continuumio/miniconda3

# Use bash instead of sh for RUN commands
# This is required for sourcing conda properly
SHELL ["/bin/bash", "-c"]

# Copy project files into container
COPY . /app

# Set working directory inside container
WORKDIR /app

# Install required system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends libltdl7 && \
    rm -rf /var/lib/apt/lists/*

# Create the main conda environment from the environment file
RUN conda env create -f natrix2.yaml && \
    conda clean -afy

# Make docker pipeline script executable
RUN chmod +x docker_pipeline.sh

# Pre-build all Snakemake conda environments during image build
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate natrix2 && \
    \
    # Create temporary dummy directory and file required for setup
    mkdir -p docker_env_build && \
    touch docker_env_build.csv && \
    \
    # Generate the dataframe (this also creates required input files like units.tsv)
    python create_dataframe.py docker_setup/env_build.yaml && \
    \
    # Create all required Snakemake conda environments only
    snakemake --configfile docker_setup/env_build.yaml --cores 1 --use-conda --conda-create-envs-only && \
    \
    # Clean up temporary files and folders
    rm -rf docker_env_build && \
    rm -f docker_env_build.csv units.tsv && \
    rm -rf results_Illumina_Swarm_16S_Prok

# Run the pipeline script with the specified project name
CMD ["sh","-c", "./docker_pipeline.sh \"$PROJECT_NAME\""]
