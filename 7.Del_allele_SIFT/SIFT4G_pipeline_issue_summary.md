# SIFT4G Annotation Pipeline — Issue Summary

## Context

Species: *Begonia masoniana* (Bmas)  
Pipeline: Liftover VCF → build snpEff database → build SIFT4G predictions → annotate variants  
Cluster: SLURM-based HPC

Reference genome has two types of sequences:
- 15 large chromosomal scaffolds: `scaffold1` – `scaffold15`
- Many smaller contigs: `Contig1453`, `Contig2353`, etc.

All variants in the VCF are on `scaffold1`–`scaffold15` only. The snpEff database (`sequences.fa`) also only contains the 15 scaffolds.

---

## The Problem

**SIFT4G_Annotator.jar requires `.regions` files named after chromosomes/scaffolds** (e.g. `scaffold1.regions`). These are produced by `sift4g` (the scorer) when the query protein FASTA headers match the scaffold names.

The GFF uses transcript IDs as feature identifiers (e.g. `ID=Bma000001.1`), and `gffread` names extracted protein sequences after those transcript IDs. So the protein FASTA has headers like `>Bma000001.1` rather than `>scaffold1`.

When `sift4g` was run with transcript ID headers, it produced **22,861 `.SIFTprediction` files** named after transcripts (e.g. `Bma000001.1.SIFTprediction`). These are a valid output format but are **not what SIFT4G_Annotator reads**. The annotator looks for `scaffold1.regions`, finds nothing, and skips all variants.

---

## What Was Tried

1. **Script 2.5** was written to rename protein FASTA headers from transcript IDs to scaffold names using a transcript → scaffold mapping derived from the GFF (column 1 of mRNA features).

2. This correctly maps e.g. `>Bma000001.1` → `>Contig1453` and `>Bma030153.1` → `>scaffold9`.

3. However, there are **22,731 proteins on scaffold sequences** distributed across the 15 scaffolds, with hundreds to thousands of transcripts per scaffold (e.g. scaffold1 has 2,733). There are also 130 proteins on Contig sequences, which are irrelevant since no variants map there.

4. Using scaffold names as protein headers means many sequences share the same header (e.g. 2,733 sequences all named `>scaffold1`). The naive approach of keeping only the longest per scaffold was rejected because it discards most of the protein information and would produce poor SIFT scores.

---

## The Core Issue — Unresolved

The fundamental tension is:

- **SIFT4G_Annotator needs one `.regions` file per scaffold**, named after the scaffold.
- **`sift4g` produces one output file per query sequence**, named after the protein header.
- **There are thousands of proteins per scaffold**, and we want SIFT scores informed by all of them, not just one representative.

The question is: **what is the correct way to provide multiple proteins per scaffold to SIFT4G while still producing scaffold-named `.regions` files that SIFT4G_Annotator can find?**

---

## Possible Avenues to Investigate

1. **Concatenate all proteins per scaffold into a single sequence** before running `sift4g`. This is biologically incorrect but may be what the tool expects for multi-gene scaffolds.

2. **Check SIFT4G documentation or GitHub issues** for how multi-transcript scaffolds are handled. The tool may have a flag or expected input format for this case.
   - SIFT4G GitHub: https://github.com/pauline-ng/SIFT4G_Annotator

3. **Check whether `.SIFTprediction` files can be renamed/symlinked to `.regions`** — it is possible the two formats are identical or compatible and only the extension differs.
   ```bash
   head -5 ~/software/snpEff/data/Bmas/Bma000001.1.SIFTprediction
   ```

4. **Use a different annotation approach** — run snpEff first to get transcript-level ANN fields, then use the transcript IDs to look up the `.SIFTprediction` files directly, bypassing SIFT4G_Annotator entirely.

---

## Current State of Files

| File | Location | Status |
|------|----------|--------|
| Lifted + sorted VCF | `$OUTPUT_DIR/all_renamed.genome_coords.sorted.vcf.gz` | ✓ Complete |
| snpEff database | `~/software/snpEff/data/Bmas/snpEffectPredictor.bin` | ✓ Complete |
| Protein FASTA (transcript headers) | `$OUTPUT_DIR/Bmas_proteins.fa` | ✓ Complete (22,861 sequences) |
| Protein FASTA (scaffold headers) | `$OUTPUT_DIR/Bmas_proteins_renamed.fa` | ✓ Produced by script 2.5 (not needed in current approach) |
| `.SIFTprediction` files | `~/software/snpEff/data/Bmas/*.SIFTprediction` | ✓ Complete (22,861 files, transcript-named) |
| `.regions` files | `~/software/snpEff/data/Bmas/*.regions` | ✗ Not produced — this is the blocker |
| SIFT-annotated VCF | `$OUTPUT_DIR/sift_results/` | ✗ Not yet produced |
| Final snpEff VCF | `$OUTPUT_DIR/all_annotated_snpeff.vcf` | ✗ Not yet produced |

---

## Scripts in the Pipeline

| Script | Purpose | Status |
|--------|---------|--------|
| `script_1_lift_VCF.sh` | Liftover VCF from bait to genome coordinates | ✓ Done |
| `script_2_build_snpeff_db.sh` | Build snpEff database for Bmas | ✓ Done |
| `script_2.5_rename_proteins.sh` | Rename protein headers (transcript ID → scaffold name) | ✓ Written, approach on hold |
| `script_3_build_regions.sh` | Run `sift4g` scorer to produce per-transcript predictions | ✓ Done (produced `.SIFTprediction`, not `.regions`) |
| `script_4_SIFT_annotate.sh` | Run SIFT4G_Annotator + snpEff annotation | ✗ Blocked — needs `.regions` files |
