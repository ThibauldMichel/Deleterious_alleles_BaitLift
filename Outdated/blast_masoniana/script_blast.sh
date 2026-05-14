#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=16G

# Load conda into this shell session
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

conda activate blast_env


# Paths
SUBJ="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
QUERY="intersected_baits.fasta"
OUTPUT="results_blast.csv"
FILT="results_filtered.csv"

# Step 1:Run blastn

blastn \
  -subject "$SUBJ" \
  -query "$QUERY" \
  -out "$OUTPUT" \
  -outfmt '6 qseqid sseqid pident length sstart send slen evalue bitscore' \
  -evalue 1e-20 \
  -max_target_seqs 50



# Step 2: Filtering for best hit per locus
# Sort by query, subject, then descending bitscore
sort -k1,1 -k2,2 -k9,9nr "$OUTPUT" |
awk '
function overlap(a_start, a_end, b_start, b_end) {
    return (a_start <= b_end && b_start <= a_end)
}

{
    q=$1; s=$2;
    start=$5; end=$6;
    if (start > end) { tmp=start; start=end; end=tmp }

    key=q"|"s;
    if (!(key in chosen)) {
        chosen[key]=start":"end;
        print;
    } else {
        split(chosen[key], ranges, " ");
        keep=1;
        for (r in ranges) {
            split(ranges[r], coords, ":");
            if (overlap(start, end, coords[1], coords[2])) {
                keep=0; break;
            }
        }
        if (keep) {
            chosen[key]=chosen[key]" "start":"end;
            print;
        }
    }
}' > "$FILT"

