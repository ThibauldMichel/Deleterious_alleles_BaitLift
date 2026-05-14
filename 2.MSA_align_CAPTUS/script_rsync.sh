#!/bin/bash
#SBATCH --job-name="rsynccapt"
#SBATCH --export=ALL
#SBATCH --mem=200G  
#SBATCH --partition=short


mkdir raw_reads

# Path to Hillebrandia reads:
rsync -azvh  /home/tmichel/projects/rbge/Begonia_genomes/Hillebrandia_illumina/hillebrandia_forward_paired.fq.gz ./raw_reads/
rsync -azvh  /home/tmichel/projects/rbge/Begonia_genomes/Hillebrandia_illumina/hillebrandia_reverse_paired.fq.gz ./raw_reads/

# Path to PNG reads:
rsync -azvh  /home/tmichel/projects/rbge/PNG_Begonia/All/* ./raw_reads/
