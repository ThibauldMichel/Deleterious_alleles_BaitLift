#!/bin/bash

#SBATCH --job-name="filtering_msa"
#SBATCH --export=ALL
#SBATCH --mem=8G
#SBATCH --partition=short

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate snakemake_env



python3 script_filtering_MSA.py
