#!/bin/bash
#SBATCH --job-name="SIFT"
#SBATCH --export=ALL
#SBATCH --mem=32G



# ===== STEP 1: Configure snpEff.config =====
# Add to snpEff/snpEff.config:
#   Bmas.genome : Bmas

# ===== STEP 2: Prepare data directory =====
mkdir -p ./snpEff/data/Bmas
mkdir -p ./snpEff/data/genomes

GFF=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff
FA=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
CDS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.cds
PEP=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.pep

cp $GFF  ./snpEff/data/Bmas/genes.gff
cp $FA       ./snpEff/data/genomes/Bmas.fa
cp $CDS        ./snpEff/data/Bmas/cds.fa       # optional but recommended
cp $PEP    ./snpEff/data/Bmas/protein.fa   # optional but recommended

# ===== STEP 3: Build the SnpEff database =====
java -Xmx8g -jar ./snpEff/snpEff.jar build -gff3 -v Bmas
# Add -noCheckCds -noCheckProtein if you don't have cds.fa/protein.fa

# ===== STEP 4: Annotate VCF with SnpEff =====
INPUT="/home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all.vcf"
SNPEFF_OUT="all_snpeff_annotated.vcf"
java -Xmx32g -jar ./snpEff/snpEff.jar -v Bmas $INPUT > $SNPEFF_OUT

# ===== STEP 5: Run SIFT4G on the SnpEff-annotated VCF =====
java -jar SIFT4G_Annotator.jar \
    -c \
    -i $SNPEFF_OUT \
    -d ./snpEff/data/Bmas \
    -r results \
    -t

# -c    Command-line mode (as opposed to GUI mode). Essential for headless environments.
# -i    Input file: path to the input VCF file (in your case, all.vcf).
# -d    Database directory: location of the SIFT4G database directory (here, CADRE.23). This folder should contain the pre-built SIFT4G database for the organism.
# -r    Results directory: directory where the annotation results will be saved (results). SIFT will generate an annotated VCF and possibly other output files here.
# -t    Use multiple threads (multithreading for faster annotation). Optional but recommended for speed.



