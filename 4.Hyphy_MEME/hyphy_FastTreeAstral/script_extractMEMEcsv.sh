#!/bin/bash

#SBATCH --job-name="csv converter"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=short





while read file;
do
	python3 script_extractMEMEcsv.py output_meme/"$file".json ;
done < mfa_files.txt
