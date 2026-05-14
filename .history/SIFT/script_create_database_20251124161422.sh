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
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.fa    data/Bmas/Bmas.fa
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.gff   data/Bmas/genes.gff
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.cds   data/Bmas/cds.fa
rsync -azvh ~/projects/rbge/tmichel/reference_genomes/Bmas.pep   data/Bmas/protein.fa

#3. Add the genome to the snpEff config file

# backup first
cp snpEff.config snpEff.config.bak

# append an entry
echo -e "\n# Begonia masoniana\nBmas.genome : Begonia_masoniana" >> snpEff.config

# 4. Build the database
cd /home/tmichel/scratch/Deleterious_alleles_PNG/SIFT/snpEff

# Run build (verbose)
###############################################
#!/bin/bash
#SBATCH --job-name="SIFT"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=200G


java -Xmx8g -jar snpEff.jar build -gff3 -v Bmas
###############################################

# 5. Troubleshooting if failure
# Database seems to be correctly built.
# Several warnings "WARNING_GENE_NOT_FOUND", as there are mRNA entries in the GFF file that are not inside the gene entries. nSnpEff therefore creates new "null" genes automatically. 