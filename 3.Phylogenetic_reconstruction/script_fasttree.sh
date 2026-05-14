#!/bin/bash
#SBATCH --job-name="gene trees PNG"
#SBATCH --export=ALL
#SBATCH --mem=8G  
#SBATCH --partition=long



#conda create -n FastTree_env -c bioconda -c conda-forge fasttree
conda activate FastTree



	mkdir -p ./05_trees/FastTree

for aln in ./04_alignments/03_trimmed/06_informed/01_coding_NUC/02_NT/*.fna; do
    base=$(basename "$aln" .fasta)
    echo "Running FastTree on $base..."
    FastTree -nt -gtr -gamma "$aln" > ./05_trees/FastTree/${base}_FastTree.tree
done







# -nt	Nucleotide sequences	Tells FastTree to treat the alignment as nucleotide data (not amino acids). It uses nucleotide substitution models.
# -gtr	General Time Reversible model	Uses the full GTR model, which estimates six substitution rate parameters (A↔C, A↔G, A↔T, C↔G, C↔T, G↔T) and empirical base frequencies — the most general reversible model commonly used in phylogenetics.
# -gamma	Gamma-distributed rate heterogeneity	Models variation in substitution rates among sites with a discrete gamma distribution. This accounts for some sites evolving faster than others — an important realism that improves branch length estimates.
