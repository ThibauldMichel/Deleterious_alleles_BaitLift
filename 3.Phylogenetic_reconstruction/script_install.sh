#!/bin/bash
#SBATCH --job-name=install_nw_prune
#SBATCH --export=ALL
#SBATCH --partition=short


# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate captus_env

#conda install -y -c conda-forge -c bioconda newick-utils


# conda install -c conda-forge biopython


conda install -c conda-forge -c bioconda amas -y

