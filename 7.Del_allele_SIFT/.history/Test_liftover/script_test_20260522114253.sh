#!/bin/bash
 
# Usage: bash compare_snippets.sh bait_locus_annotation.tsv baits.fasta genome.fasta
 
TSV="/home/tmichel/scratch/Deleterious_alleles_BaitLift/1.Baits_to_genome/bait_locus_annotation.tsv"
BAITS="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
GENOME="/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta"
WINDOW=50
 
# Skip header, parse all 18 columns explicitly
tail -n +2 "$TSV" | while IFS=$'\t' read -r scaffold t_start t_end bait_id mapq strand q_start q_end q_len \
    gff_seqid gff_source gff_type gff_start gff_end gff_score gff_strand gff_phase gff_attributes; do
 
    # Extract bait snippet (samtools faidx uses 1-based coords)
    bait_snip=$(samtools faidx "$BAITS" "${bait_id}:$((q_start+1))-$((q_start+WINDOW))" 2>/dev/null | grep -v "^>" | tr -d '\n')
 
    # Extract genome snippet
    genome_snip=$(samtools faidx "$GENOME" "${scaffold}:$((t_start+1))-$((t_start+WINDOW))" 2>/dev/null | grep -v "^>" | tr -d '\n')
 
    # Reverse complement genome snippet if minus strand
    if [ "$strand" == "-" ]; then
        genome_snip=$(echo "$genome_snip" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    fi
 
    echo "---"
    echo "Bait  : $bait_id  q_coords: ${q_start}-$((q_start+WINDOW))  strand: $strand"
    echo "Locus : ${scaffold}:${t_start}-${t_end}"
    echo "Bait    : $bait_snip"
    echo "Genome  : $genome_snip"
 
done
