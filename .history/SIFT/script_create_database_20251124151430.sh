#!/bin/bash


# Building a new database on the B. masoniana reference genome.
# https://pcingola.github.io/SnpEff/snpeff/build_db/

# Step 1: configure a new genome


# The working directory:
/home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/snpEff

# The database directory:
~/projects/rbge/tmichel/reference_genomes

# The database files:
Bmas.cds 
Bmas.fa 
Bmas.gff 
Bmas.pep

#Bmas.fa — genome FASTA
#Bmas.gff — annotation GFF (GFF3)
#Bmas.cds — CDS FASTA
#Bmas.pep — protein FASTA

# 1. Create snpEff data directory for the new genome

# Go to snpEff dir
cd /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/snpEff

# Create the data folder for the genome
mkdir -p data/Bmas

# 2. Copy the genome files to the new data folder

# Copy your files into snpEff structure and rename them to the expected names
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.fa    data/genomes/Bmas.fa
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.gff   data/Bmas/genes.gff
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.cds   data/Bmas/cds.fa
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.pep   data/Bmas/protein.fa
