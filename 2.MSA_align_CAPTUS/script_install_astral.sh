#!/bin/bash
#SBATCH --job-name="install ASTRAL"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=medium


source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate FastTree_env

conda install -c bioconda astral --yes
