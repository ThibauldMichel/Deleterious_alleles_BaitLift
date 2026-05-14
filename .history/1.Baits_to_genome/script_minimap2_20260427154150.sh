#!/bin/bash
#SBATCH --job-name="minimap2"
#SBATCH --export=ALL
#SBATCH --mem=16G  
#SBATCH --partition=medium

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate bcftools_env

BAITS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Hannah_Begonia_baits_edited.fasta
REFERENCE=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
GFF=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff

# ── PAF output: query + target coordinates ─────────────────────────────────────
# PAF (Pairwise Alignment Format) natively reports both query (bait) and target
# (genome) coordinates in the same line. No -a flag — minimap2 outputs PAF by
# default without it.
minimap2 -x splice \
         --cs \
         -C5 \
         -u b \
         $REFERENCE \
         $BAITS \
         > bait_hits.paf

# -x splice  Splice-aware preset: allows large gaps corresponding to introns between
#            exon blocks. Appropriate here because baits are CDS-derived (B. luzhaiensis
#            transcriptome) and will span introns when mapped to the genomic reference.
# --cs       Adds a cs tag encoding base-level alignment differences (mismatches,
#            insertions, deletions). Useful for fine-scale mismatch info and downstream
#            parsing of sequence variation.
# -C5        Reduces the penalty for non-canonical splice sites, important because
#            splice sites may differ slightly between B. luzhaiensis and B. masoniana.
# -u b       Search for GT-AG splice motifs on both strands rather than inferring
#            strand from sequence alone — sensible for cross-species mapping where
#            strand signal may be weaker.

# PAF output columns:
# col 1   query name       bait contig name
# col 2   query length     total bait sequence length
# col 3   query start      0-based position on the bait where the match begins
# col 4   query end        position on the bait where the match ends
# col 5   strand           + or - relative to the reference
# col 6   target name      scaffold / chromosome name
# col 7   target length    total scaffold length
# col 8   target start     genome coordinate (0-based)
# col 9   target end       genome coordinate
# col 10  residue matches  number of matching bases
# col 11  block length     total alignment block length including gaps
# col 12  mapq             mapping quality (60 = unique best hit)

# ── Filter: keep only primary alignments (mapq == 60) ─────────────────────────
# PAF does not use SAM flags, but mapq = 0 indicates secondary or ambiguous hits.
# Filtering on mapq >= 60 retains only unique, high-confidence primary alignments,
# equivalent to -F 2308 in the SAM/BAM pipeline.
awk '$12 >= 60' bait_hits.paf > bait_hits_primary.paf

# ── Filter: keep only baits with >= 70% query coverage ────────────────────────
# Baits where less than 70% of the bait sequence aligned indicate only partial
# anchoring to the genome. These are excluded from downstream GERP/HyPhy/SIFT
# analyses as the genomic locus may not be reliably identified.

# Sort PAF by baits and query coordinates
sort -k1,1 -k3,3n bait_hits_primary.paf > bait_hits_primary.sorted.paf

# Total coverage per bait
awk '
BEGIN { OFS="\t" }

function flush() {
    if (bait != "") {
        cov_pct = total_cov / bait_len * 100
        print bait, bait_len, total_cov, sprintf("%.1f", cov_pct)
    }
}

{
    if ($1 != bait) {
        flush()
        bait = $1
        bait_len = $2
        total_cov = 0

        start = $3
        end   = $4
    } else {
        if ($3 <= end) {
            if ($4 > end) end = $4
        } else {
            total_cov += (end - start)
            start = $3
            end   = $4
        }
    }
}

END {
    total_cov += (end - start)
    flush()
}
' bait_hits_primary.sorted.paf > bait_coverage.tsv

# Filter bait by 70% total coverage
awk '$4 >= 70' bait_coverage.tsv > bait_70pct_ids.txt

# Keep only those baits in original PAF
grep -Ff <(cut -f1 bait_70pct_ids.txt) bait_hits_primary.paf > bait_hits_70pct_merged.paf



# ── Convert PAF to extended TSV with both query and target coordinates ─────────
# PAF columns used:
#   $1  = query name (bait contig)
#   $2  = query length
#   $3  = query start (0-based)
#   $4  = query end
#   $5  = strand
#   $6  = target name (scaffold)
#   $8  = target start (0-based)
#   $9  = target end
#   $10 = residue matches
#   $11 = alignment block length
#   $12 = mapq
# Query coverage and identity are computed on the fly.

awk 'BEGIN{
    OFS="\t"
    print "bait_name", "bait_length", "bait_start", "bait_end", \
          "strand", "scaffold", "target_start", "target_end", \
          "matches", "block_length", "mapq", "query_cov_pct", "identity_pct"
}
{
    qcov = ($4 - $3) / $2 * 100
    ident = $10 / $11 * 100
    print $1, $2, $3, $4, \
          $5, $6, $8, $9, \
          $10, $11, $12, \
          sprintf("%.1f", qcov), sprintf("%.1f", ident)
}' bait_hits_70pct_merged.paf > bait_hits_70pct_extended.tsv


