#!/bin/bash
#SBATCH --job-name=rename_msa
#SBATCH --export=ALL
#SBATCH --mem=8G
#SBATCH --partition=short


ALIGN="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/04_alignments/03_trimmed/06_informed/01_coding_NUC/03_genes"

mkdir ./MSA_renamed


for aln in "$ALIGN"/*.fna; do
 gene=$(basename "$aln" .fna) ;
    sed 's/ .*//' "$aln" > ./MSA_renamed/${gene}.fna ;
done
