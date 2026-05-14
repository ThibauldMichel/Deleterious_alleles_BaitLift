#!/bin/bash
#SBATCH --job-name="ASTRAL species tree"
#SBATCH --export=ALL
#SBATCH --mem=64G  
#SBATCH --partition=long

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate FastTree_env

# Make a folder for species tree
mkdir -p ./05_phylogeny_FastTreeAstral/ASTRAL

# concatenate (assumes one tree per file)
cat ./05_phylogeny_FastTreeAstral/genes_tree_FastTree/*.tree > ./05_phylogeny_FastTreeAstral/ASTRAL/gene_trees_all.tre

# quick sanity check: count trees
echo "Number of gene trees:" $(wc -l < ./05_phylogeny_FastTreeAstral/ASTRAL/gene_trees_all.tre)

# draw the Astral tree:
java -Xmx32G -jar  /home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/software/ASTRAL/Astral/astral.5.7.8.jar \
	 -i ./05_phylogeny_FastTreeAstral/ASTRAL/gene_trees_all.tre \
      -o ./05_phylogeny_FastTreeAstral/ASTRAL/species_tree_astral.tre \
      -t 3 \
     > ./05_phylogeny_FastTreeAstral/ASTRAL/astral_run.log 2>&1



# java -Xmx32G	Allocates 32 GB of memory for the Java Virtual Machine running ASTRAL. Adjust this based on your available resources.
# -jar /path/to/astral.5.7.8.jar	Specifies the ASTRAL executable JAR file to run.
# -i ./gene_trees_all.tre	Input file containing all individual gene trees (one per line or concatenated). ASTRAL will use these to infer the species tree.
# -o ./species_tree_astral.tre	Output file where the inferred species tree will be written.
# -t 3	Requests that ASTRAL output additional information — this is the “output option flag”, where:
	# -t 1 → only species tree (default)
	# -t 2 → species tree + local posterior probabilities (LPPs)
	# -t 3 → species tree + branch quartet scores + local posterior probabilities (the most detailed output)
# > astral_run.log 2>&1	Redirects both standard output (stdout) and standard error (stderr) into the file astral_run.log, so the run’s full log is saved.
