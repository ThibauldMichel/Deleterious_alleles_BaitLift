#!/bin/bash

#SBATCH --job-name="mergecsv"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=short


python3 script_merge_meme_csv.py


