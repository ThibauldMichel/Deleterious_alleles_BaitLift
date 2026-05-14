#!/bin/bash
#SBATCH --job-name="captus PNG"
#SBATCH --export=ALL
#SBATCH --mem=500G  
#SBATCH --partition=himem

#conda create -n captus -c bioconda -c conda-forge captus iqtree
#conda activate captus

# 0. Check installation
#captus -h

# 1. Cleaning Reads
#captus clean -r ./raw_reads

# 2. De Novo Assembly
#captus assemble -r 01_clean_reads --max_contig_gc 70

# 3. Extracting Target Sequences
#captus extract -a 02_assemblies -n /home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta

# 4. Multiple Sequence Alignment
#captus align -e 03_extractions

# 5. Phylogenetic inference
#mkdir 05_phylogeny && cd 05_phylogeny
#iqtree -p ../04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT -pre concat -T AUTO --single-output



# Step 1: Run ModelFinder Only

 iqtree -p ./04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT -pre concat_modeltest -m MFP -T AUTO -keep-ident 

# -p for partitioned file.
# -m MFP activates ModelFinder Plus to test many models.
 
# It will output a file like: concat_modeltest.best_scheme.nex and a log in concat_modeltest.log.

# It performs a partitioned phylogenetic analysis on nucleotide alignments stored in the specified directory, automatically chooses the best model for each partition, and constructs a maximum-likelihood tree. It keeps original sequence names and optimizes performance with automatic multi-threading.
# In a phylogenetic context, partitions refer to subsets of your alignment that are treated independently when estimating the model of evolution.
# Each partition:
# -Can be a different gene or genomic region,
# -Can evolve under a different evolutionary model,
# -May have different rates or patterns of mutation.


#??????????????????????????????????????????????????????????????????
# How to enable resumable runs?

# The command you posted already generates a checkpoint file automatically:

# concat_modeltest.ckp.gz

# This file stores the state of the run at the last checkpoint.
# To resume a run after interruption:

# Just rerun the exact same command, and IQ-TREE will detect the checkpoint file and resume automatically:
# iqtree -p ./04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT -pre concat_modeltest -m MFP -T AUTO -keep-ident
# Q-TREE will see that concat_modeltest.ckp.gz exists and resume instead of starting over.

#??????????????????????????????????????????????????????????????????






# Step 2: Use This Model in Final Tree Search

#iqtree -p ./04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT -pre concat_fixedmodel -m GTR+F+R5 -T AUTO --single-output -keep-ident

iqtree2 \
    -p ./04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT \
    -spp concat_modeltest.best_scheme.nex \
    -pre concat_finaltree \
    -T AUTO \
    --single-output \
    -keep-ident \
    -bb 1000 \
    -alrt 1000

# -T AUTO
# Tells IQ-TREE to automatically determine the number of CPU threads to use (usually all available cores).

# --single-output
# Consolidates output into a single directory (instead of creating multiple subfolders for each partition or model).

# -keep-ident
# Prevents IQ-TREE from modifying sequence names (e.g., by shortening them or replacing special characters).

# -bb 1000: runs 1000 ultrafast bootstrap replicates.

# -alrt 1000: runs 1000 SH-aLRT replicates, a fast approximation of branch support.

# -spp: tells IQ-TREE to use your partition scheme with specific models per partition (which is what *.best_scheme.nex contains).

