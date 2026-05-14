#!/bin/bash
#SBATCH --job-name=iqtree_genetrees
#SBATCH --export=ALL
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --array=0-1049%20        # 1050 genes total, run max 50 in parallel
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err
#SBATCH --cpus-per-task=4


# --- Activate environment ---
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate captus_env

# --- Define directories ---
MSA="/home/tmichel/scratch/Deleterious_alleles_PNG/CAPTUS/CAPTUS-PNG/04_alignments/03_trimmed/06_informed/01_coding_NUC/03_genes"
OUTDIR="./05_iqtree_genes"
mkdir -p "$OUTDIR" logs

# --- Select the current alignment ---
FILES=("$MSA"/*.fna)
aln="${FILES[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$aln" .fna)

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID starting on $base"

# --- Skip if output tree already exists (safety checkpoint) ---
if [[ -f "$OUTDIR/${base}.treefile" ]]; then
    echo "Tree for $base already exists. Skipping."
    exit 0
fi

# --- Run IQ-TREE2 ---
iqtree3 \
  -s "$aln" \
  -m MFP \
  -T $SLURM_CPUS_PER_TASK \
  -pre "$OUTDIR/$base" \
  -keep-ident \
  -redo

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID finished $base"



# -s "$aln"	Sequence alignment file	Specifies the input multiple sequence alignment (MSA) file for one gene. It can be in FASTA, PHYLIP, or NEXUS format. IQ-TREE reads this file and builds an ML tree from it.
# -m MFP	ModelFinder Plus	Tells IQ-TREE to automatically test a wide range of nucleotide substitution models (e.g., JC, HKY, TN, GTR, with +I, +G, +I+G) and select the best-fitting model for that alignment using AICc or BIC. This gives you the most appropriate model per gene without having to specify it manually.
# -T $SLURM_CPUS_PER_TASK	Number of CPU threads	Informs IQ-TREE how many CPU cores to use for this job. The variable $SLURM_CPUS_PER_TASK comes from your SLURM submission (#SBATCH --cpus-per-task=4), ensuring each array job uses 4 threads. You can also use -T AUTO to let IQ-TREE detect cores automatically.
# -pre "$OUTDIR/$base"	Output file prefix	All output files (e.g., .log, .treefile, .iqtree, .model.gz, etc.) will be written to $OUTDIR and prefixed with the gene name ($base). Example: ./05_trees/IQ/gene123.treefile. This keeps each gene’s results organized.
# -keep-ident	Preserve sequence identifiers	Prevents IQ-TREE from shortening or modifying long taxon names in the alignment. Without this flag, IQ-TREE may truncate long names for internal handling (e.g., B.S.Pet.ELAE86.218 → B_S_Pet_ELAE86_218). Keeping identifiers intact ensures consistency with other tools like HyPhy or GERP.
# -redo	Overwrite existing output	Tells IQ-TREE to overwrite previous results if a run with the same -pre prefix already exists. Without this, IQ-TREE would refuse to run again unless you delete the old files manually. Very handy for reruns and job arrays.
