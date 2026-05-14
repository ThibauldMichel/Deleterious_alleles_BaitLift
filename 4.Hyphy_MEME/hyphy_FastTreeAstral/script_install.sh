#!/bin/bash

#SBATCH --job-name="install biopython"
#SBATCH --export=ALL
#SBATCH --mem=64G
#SBATCH --partition=short

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate hyphy_env



#conda create -n hyphy_env -c bioconda hyphy
#conda create -n newick_env bioconda::newick_utils
#conda activate newick_env
#conda install bioconda::newick_utils


conda install -c conda-forge biopython
