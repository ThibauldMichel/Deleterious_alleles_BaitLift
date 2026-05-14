# BaitSift

**A pipeline for detecting deleterious alleles from Hyb-Seq bait capture data**

BaitSift bridges targeted sequence capture (Hyb-Seq) assemblies to a reference genome, enabling exon-aware coordinate liftover, multi-sequence alignment, phylogenetic reconstruction, and functional annotation of deleterious coding variants. The pipeline was developed to characterise the load of deleterious alleles across *Begonia* samples using a bait set derived from the *B. luzhaiensis* transcriptome.

---

## Pipeline overview

```
Hyb-Seq assemblies (CAPTUS)
        │
        ▼
1. Baits to genome          — Align bait sequences to the B. masoniana reference genome
        │                      using minimap2 (-ax splice); resolve exon block structure;
        │                      produce a coordinate mapping table (bait space → genome space)
        ▼
2. MSA alignment (CAPTUS)   — Multiple sequence alignment of per-locus assemblies across
        │                      all samples
        ▼
3. Phylogenetic              — Reconstruct gene trees and species-level phylogeny from
   reconstruction            aligned loci
        │
        ▼
4. HyPhy MEME               — Test for episodic positive / diversifying selection at
        │                      individual codon sites (Mixed Effects Model of Evolution)
        ▼
5. GERP conservation        — Compute per-site conservation scores across a multi-species
   scores                    alignment; identify evolutionarily constrained positions
        │
        ▼
6. Variant calling &        — Call variants from Hyb-Seq assemblies in bait-relative
   annotations               coordinates; lift over to genome coordinates using the
        │                      exon-aware mapping table from Step 1
        ▼
7. Deleterious allele       — Annotate lifted variants with SIFT4G predictions;
   annotation (SIFT)         integrate SIFT scores with GERP and HyPhy MEME results
                               to identify high-confidence deleterious alleles
```

---

## Key design decisions

**Exon-aware liftover.** Bait sequences are transcriptome-derived CDS contigs. When mapped to the genomic reference with `minimap2 -ax splice`, each bait aligns across multiple exon blocks separated by introns. Variant positions in bait-relative coordinates therefore cannot be converted to genome coordinates with a simple linear offset — the liftover must resolve which exon a variant falls in and apply the correct strand-specific offset. This is handled by `liftover_bait_to_genome.py`.

**Functional annotation via SIFT4G + snpEff.** SIFT4G is used to predict the functional impact of coding variants. Annotations are generated via the snpEff pipeline (rather than SIFT4G native tools) to obtain richer outputs including transcript consequences, gene IDs, and splice effects.

**Multi-evidence prioritisation.** High-confidence deleterious variants are defined by convergent evidence across three sources: a SIFT score < 0.05 (predicted deleterious), a GERP RS score above a conservation threshold (evolutionarily constrained site), and a HyPhy MEME result consistent with purifying rather than positive selection.

---

## Reference genome & bait set

| Resource           | Details                                                                                                                    |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------- |
| Reference genome   | *Begonia masoniana* (`Bmas.fa` / `Bmas.gff`)                                                                         |
| Bait source        | *B. luzhaiensis* transcriptome (Tseng et al., 2017), filtered against *Cucumis sativus* annotation (Yang et al., 2012) |
| Bait set additions | Angiosperms353 matches (Johnson et al., 2019); developmental genes; light-response pathway genes                           |
| Mapped baits       | 1,218 / 1,250 baits successfully mapped to*B. masoniana* (21 unmapped)                                                   |

---

## Directory structure

```
BaitSift/
├── 1.Baits_to_genome/
│   ├── script_minimap2.sh
│   ├── script_check_multiple_location.sh
│   ├── script_finding_missing_baits.sh
│   └── bait_locus_annotation.tsv        ← coordinate mapping table
├── 2.MSA_align_CAPTUS/
├── 3.Phylogenetic_reconstruction/
├── 4.Hyphy_MEME/
├── 5.GERP_cons_scores/
├── 6.Variant_calling_annotations/
│   └── all_renamed.vcf.gz               ← variants in bait-relative coords
└── 7.Del_allele_SIFT/
    ├── liftover_bait_to_genome.py
    ├── sift_annotation_workflow.sh
    └── LIFTOVER_AND_SIFT_GUIDE.md
```

---

## Dependencies

- `minimap2`
- `samtools`
- `bedtools`
- `bcftools`
- `CAPTUS`
- `HyPhy`
- `SIFT4G` (SIFT4G_Create_Genomic_DB, SIFT4G_Annotator)
- `snpEff`
- Python 3.6+ (`pandas`)
