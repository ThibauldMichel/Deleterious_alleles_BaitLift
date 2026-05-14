#!/bin/bash
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate snakemake_env


python3 script_sync_msa_tree.py
