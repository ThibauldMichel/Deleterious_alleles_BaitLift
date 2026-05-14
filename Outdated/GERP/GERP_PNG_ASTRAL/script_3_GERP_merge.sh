#!/bin/bash
#SBATCH --job-name="GERP_merge"
#SBATCH --output=merge_%j.out
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --partition=short



cat mfa_files.txt | while read file; do
    awk -v name="$file" '{print name, $0}' OFS='\t' "mfa_dir/${file}_edited.mfa.rates.elems"
done > all.mfa.rates.elems

