#!/bin/bash
#SBATCH --job-name="install"
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --export=ALL


# Load conda into this shell session
#source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

#conda activate snakemake_env



#mamba install -y -c bioconda newick_utils
conda create -n captus_env -c conda-forge -c bioconda iqtree --y



