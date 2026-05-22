#!/bin/bash
 
# Usage: bash compare_snippets.sh bait_locus_annotation.tsv baits.fasta genome.fasta
 
TSV="/home/tmichel/scratch/Deleterious_alleles_BaitLift/1.Baits_to_genome/bait_locus_annotation.tsv"
BAITS="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
GENOME="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
WINDOW=50
 
# Skip header
tail -n +2 "$TSV" | while IFS=$'\t' read -r scaffold t_start t_end bait_id mapq strand q_start q_end q_len rest; do
 
    # Extract bait snippet (samtools faidx uses 1-based coords)
    bait_snip=$(samtools faidx "$BAITS" "${bait_id}:$((q_start+1))-$((q_start+WINDOW))" | grep -v "^>" | tr -d '\n')
 
    # Extract genome snippet
    genome_snip=$(samtools faidx "$GENOME" "${scaffold}:$((t_start+1))-$((t_start+WINDOW))" | grep -v "^>" | tr -d '\n')
 
    # Reverse complement genome snippet if minus strand
    if [ "$strand" == "-" ]; then
        genome_snip=$(echo "$genome_snip" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    fi
 
    echo "---"
    echo "Bait  : $bait_id  coords: ${q_start}-$((q_start+WINDOW))  strand: $strand"
    echo "Locus : ${scaffold}:${t_start}-${t_end}"
    echo "Bait    : $bait_snip"
    echo "Genome  : $genome_snip"
 
done
