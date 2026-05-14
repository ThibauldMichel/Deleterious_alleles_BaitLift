#!/bin/bash
#SBATCH --job-name="SIFT"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=200G


java -Xmx8g -jar snpEff.jar build -gff3 -v Bmas
