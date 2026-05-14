# VCF Liftover and SIFT4G Annotation Guide

## Overview

Your VCF file contains variants called in **bait-relative coordinates** (from CAPTUS assembly), but SIFT4G requires **genome coordinates** (to map variants to annotated CDS features in your GFF).

The mapping table from Stage 1 contains all the information needed to convert positions between these coordinate systems.

---

## Key insight: Exon-aware position mapping

Your baits are **transcriptome-derived CDS sequences** that span **multiple exons** when mapped to the genome. For example:

```
Bait: Hannah_Begonia_baits-ACmerged_contig_1002 (811 bp)
  └─ Maps to scaffold2:49251268–49439930 (+)
      └─ Exon 1: scaffold2:49252054–49252134 (query 0–80)
      └─ Exon 2: scaffold2:49252905–49253100 (query 80–275)
      └─ Exon 3: scaffold2:49253357–49253472 (query 275–390)
      └─ Exon 4: scaffold2:49256999–49257301 (query 390–692)
      └─ ... (more exons)
```

A variant at **position 150 in the bait** does NOT map to a simple linear offset. Instead:

1. Determine which exon contains position 150 (exon 2 in this example)
2. Calculate offset within that exon (150 - 80 = 70 bp)
3. Map to the corresponding genome position (49252905 + 70 = 49252975)

---

## Workflow

### Step 1: Prepare inputs

Ensure you have:

1. **VCF file** in bait-relative coordinates:

   ```
   /home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all_renamed.vcf.gz
   ```
2. **Mapping table** from Stage 1 (your coordinate map):

   ```
   Chr | t_start | t_end | bait_id | mapq | strand | q_start | q_end | q_len | gff_seqid | gff_source | gff_type | gff_start | gff_end | ...
   ```

   This should be a TSV with headers. If you created it as a bedtools output, ensure headers are added.
3. **Reference genome** and **GFF**:

   ```
   /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
   /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff
   ```

### Step 2: Run the liftover

```bash
VCF=/home/tmichel/scratch/Deleterious_alleles_PNG_baits/6.Variant_calling_annotations/all_renamed.vcf.gz
TSV=/home/tmichel/scratch/Deleterious_alleles_PNG_baits/1.Baits_to_genome/bait_locus_annotation.tsv
FAS=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
GFF=/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff
OUT=./
# Run the full workflow
bash sift_annotation_workflow.sh $VCF $TSV $FAS $GFF $OUT
```

This will:

1. Lift variant positions from bait space to genome space
2. Validate the lifted VCF
3. Build a SIFT4G database from your reference and GFF
4. Annotate all variants with SIFT scores

### Step 3: Inspect results

Key output files:

```
output_directory/
├── all_renamed.genome_coords.sorted.vcf.gz
│   └─ Lifted VCF in genome coordinates (indexed)
├── sift_results/
│   ├── SIFT_results.txt
│   │   └─ Summary of SIFT4G run
│   ├── *_exonicSIFT.vcf
│   │   └─ Annotated VCF with SIFT predictions
│   └── *_log.txt
│       └─ Detailed run log
└── sift4g_bmas/
    └─ SIFT4G database (can delete after annotation if disk space is tight)
```

### Step 4: Parse SIFT results

SIFT4G will annotate the VCF with these INFO fields:

```
SIFT=deleterious(0.01) | tolerated(0.15) ...
```

Extract deleterious variants:

```bash
bcftools filter -i 'INFO/SIFT ~ "deleterious"' \
    output_directory/sift_results/*_exonicSIFT.vcf \
    -O z -o deleterious_variants.vcf.gz
```

---

## Understanding the liftover logic

### Forward strand (+)

Position in bait → Position in genome (straightforward offset)

```
Bait query:    [====VARIANT====]
               0                 811 bp

Exon 2 maps to genome:
  Query: 80–275 (in bait)
  Genome: 49252905–49253100 (in scaffold2)

Variant at bait position 150:
  Offset in exon: 150 - 80 = 70
  Genome position: 49252905 + 70 = 49252975 ✓
```

### Reverse strand (-)

Position in bait → Position in genome (flipped direction)

```
Bait query (forward):    [====VARIANT====]
                         0                 811 bp

Exon maps to genome (reverse strand):
  Query: 80–275 (in bait, still in forward orientation)
  Genome: 49252905–49253100 (but read right-to-left in genome)

Variant at bait position 150:
  Offset in exon: 150 - 80 = 70
  Genome position (reversed): 49253100 - 70 - 1 = 49253029 ✓
```

The script handles this automatically via the `strand` field.

---

## Troubleshooting

### Problem 1: "Position not found in any exon"

This means a variant in your VCF does not align to any exon block in the mapping table.

**Causes:**

- Variant position is outside the mapped bait region
- Bait sequence was not successfully mapped to the genome (not in mapping table)
- Bait ID in VCF doesn't match the mapping table

**Solution:**
Check the failed count in the script output. If it's high:

```bash
# Re-run with more debugging
# The script will skip unmappable variants silently

# Instead, save them and investigate
```

### Problem 2: VCF validation fails after liftover

```bash
bcftools view -H output_directory/all_renamed.genome_coords.sorted.vcf.gz | head

# Check that chromosome names match your reference:
bcftools view output_directory/all_renamed.genome_coords.sorted.vcf.gz | \
    bcftools query -f '%CHROM\n' | sort -u | head
```

Compare against:

```bash
grep "^[^#]" /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff | \
    awk '{print $1}' | sort -u
```

If they don't match, check the `chr` column in your mapping table.

### Problem 3: SIFT4G database build fails

Common issues:

1. **GFF format problem**: SIFT4G is picky about GFF format. Ensure:

   ```bash
   # Check GFF has CDS features
   grep "CDS" /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff | head

   # Check column count (must be 9)
   head -1 /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff | \
       awk -F'\t' '{print NF}'  # Should output 9
   ```
2. **Sequence ID mismatch**: Ensure sequence IDs in GFF match FASTA headers:

   ```bash
   # GFF sequence IDs
   awk -F'\t' '{print $1}' /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff | \
       sort -u | head

   # FASTA headers
   grep "^>" /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa | head
   ```
3. **Indexed FASTA**: SIFT4G needs a `.fai` index:

   ```bash
   samtools faidx /home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa
   ```

---

## Integration with Stage 8 (Candidate prioritisation)

Once SIFT annotation is complete, integrate results with GERP and HyPhy MEME:

```python
import pandas as pd

# Load SIFT results
sift = pd.read_csv('sift_results/SIFT_results.txt', sep='\t', low_memory=False)

# Load GERP scores (from Stage 5)
gerp = pd.read_csv('all.mfa.rates.elems', sep='\t', usecols=['index', 'RS'])

# Load HyPhy MEME results (from Stage 4)
meme = pd.read_csv('csv_extracted_MEME/all_genes.csv', sep='\t')

# Merge on genomic position
candidates = sift.merge(gerp, left_on='position', right_on='index')
candidates = candidates.merge(meme, on='codon')

# Filter for high-confidence deleterious variants:
# SIFT < 0.05, GERP RS > threshold, omega1 < 1
high_conf = candidates[
    (candidates['SIFT_score'] < 0.05) &
    (candidates['RS'] > 3) &
    (candidates['omega1'] < 1)
]

print(f"High-confidence candidates: {len(high_conf)}")
high_conf.to_csv('deleterious_candidates.csv', index=False)
```


---

## The SIFT/snpEff issue

What tool to choose?

Two pipelines have been built: One building database with:
SIFT4G_Create_Genomic_DB (build species database) -> SIFT4G_Annotator.jar

---



## Citation & methods text

If you use this approach, consider adding to your methods:

> "Variants were called in bait-relative coordinates from Hyb-Seq assemblies (CAPTUS) and lifted to reference genome coordinates using a position-aware mapping table derived from minimap2 alignments of bait sequences to the *B. masoniana* reference genome. SIFT4G was then used to predict the functional impact of coding variants using the reference GFF annotation."

---

## Files provided

1. **`liftover_bait_to_genome.py`** – Core liftover logic
2. **`sift_annotation_workflow.sh`** – End-to-end workflow script
3. **`LIFTOVER_AND_SIFT_GUIDE.md`** – This guide

All scripts require:

- Python 3.6+
- bcftools
- SIFT4G tools (SIFT4G_Create_Genomic_DB, SIFT4G_Annotator)
- samtools (for FASTA indexing)
