
# Bridge baits set to the B.masoniana genome

## 1. The alignment of the baits to the genome. Strugglings to avoid overlapping coordinates and choosing the right options for the baits which are concatenation of exons.

We are running the script_minimap2.sh to bridge between the baits set and the B.masoniana genome.

To make the script addequate to the minimap alignment, we have inquired in the bait origin:

`We have produced a bait set that works across Begonia and includes sequences for genes likely to be of interest in the genus. We started with likely single-copy genes from a transcriptome produced from the Asian species Begonia luzhaiensis (Tseng et al., 2017). We were concerned that the baits ought not to be anonymous, because for many downstream analyses it is important to understand the function of the genes used. At the time we were designing the baits, the closest related species with a well-annotated genome available was cucumber (Cucumis sativus L.), so we limited target sequences to those annotated in the cucumber genome assembly (Yang et al., 2012) (see Figure 1). We wished to link the bait set to work already done and to maximise our ability for use in functional studies. We added sequences that matched markers used to generate the first Begonia genetic map (Brennan et al., 2012), along with genes linked to quantitative trait loci from an analysis of species-level variation (Twyford et al., 2014) and a DESeq analysis (Emelianova et al., 2021) and sequences from transcription factors differentially expressed between Begonia conchifolia and B. plebeja (Emelianova, 2017). Most Begonia are shade-adapted, and we wished to allow sequence analysis of key genes in the pathways of light perception and response. We added the Begonia luzhaiensis orthologues of genes associated with light responses, particularly shade tolerance. Based on performance in the Begonia sect. Coelocentrum capture (see Table 1), we refined the bait set, removing baits that captured many paralogues or that captured poorly, and replacing them with sequences of developmental genes and matches to the Angiosperms353 baits. We identified 13 target sequences with high paralogy and 83 that failed to capture. These sequences were removed from the set and replaced with sequences of several more developmental gene sequences matching the Angiosperms353 bait set (Johnson et al., 2019). This set was used for the three further captures that generated the data analysed here.`


*This description clarifies everything. Your baits are **transcriptome-derived CDS sequences** (from Begonia luzhaiensis transcriptome, filtered against the cucumber genome annotation, plus Angiosperms353 matches). This confirms that `-ax splice` is the correct preset — the baits are continuous exonic sequences that will span introns when mapped to the Begonia masoniana genome.*

*There is one additional nuance though. Your baits come from a **divergent transcriptome** ( B. luzhaiensis , an Asian species) being mapped to what is presumably a different Begonia species genome (Bmas.fa — likely  B. masoniana ?). The genetic distance between the bait source and your reference matters for minimap2 preset tuning.*

*The `-ax splice` preset uses relatively strict parameters optimised for mapping a transcript back to its own genome. For cross-species mapping at the within-genus level (which is your situation), you want to relax the mismatch tolerance. The recommended adjustmejnt include key additions that are `-C5` which reduces the penalty for non-canonical splice sites (important because splice sites may differ slightly between B. luzhaiensis and  B. masoniana ), and `-u b` which tells minimap2 to search for GT-AG splice motifs on both strands rather than inferring strand from the sequence alone — sensible when mapping across species where strand signal may be weaker.*

*However, looking at your BED output, your mapping quality scores are all 60 (maximum mapq), which suggests the within-genus divergence is not causing any real alignment difficulty. Your `-ax splice` script as written is already working well. The duplicate-line issue described in my previous message (same contig appearing with both `+` and `-` strand) is the only genuine problem to fix, and that is purely a `bamtobed` issue, not a minimap2 preset issue.*

*So to directly answer your original question:  **`-ax splice` is correct for your data** . The baits are CDS-derived, the mapping is cross-species within the same genus, and your mapq=60 results confirm it is working. The only change needed to your script is the `samtools view -F 2308` filter and the `-split` flag on `bamtobed` before intersecting with the GFF.*


The results of this script give favorable outputs:


*This looks much better. The duplicate strand-flip artefact is gone — every block for a given contig now has a consistent strand (`+` or `-`), which confirms the `-F 2308` filter worked correctly.*

*What you're seeing now is the  **correct biological signal** : each bait contig maps to multiple exon blocks separated by introns. For example `contig_7251` maps to 8 separate exon blocks across the region `460618–463131` on scaffold1, with gaps between blocks corresponding to introns. This is exactly what you expect from CDS-derived baits mapped to a genomic reference — the `-split` flag is correctly decomposing the spliced alignment into its constituent exon intervals.*

*This output is now ready for `bedtools intersect` with the GFF. But before that, let's check the mapping in different location.*


## 2. Check the baits are located on the same contig and not scattered across the genome

*One thing worth checking before you proceed though — you will likely have some baits that mapped to multiple genomic locations (paralogs, gene family members). It is worth quantifying this with `script_check_multiple_location.sh`. If a bait maps to several unrelated loci (not just multiple exon blocks of the same gene), that bait was flagged in your paper as a high-paralogy case and should probably be excluded before the GERP/HyPhy analysis. The `-F 2308` filter already kept only the best primary alignment per bait, so you are protected against that to some extent — but it is worth a sanity check given that you explicitly removed high-paralogy baits during bait set refinement.*

The results of script_check_multiple_location.sh are reassuring, no-nonsense locations mapped.

*The numbers here are reassuring on the paralogy front, but the exon count is worth thinking about.*

*The `1218` vs `0` result tells you that all 1,218 baits that appear in the BED have a unique name — no bait mapped to two completely separate genomic loci. The `-F 2308` primary-alignment filter did its job: each bait has one best mapping location. You have no paralogy problem to worry about.*

*The high counts (35, 33, 29...) are not multiple loci — they are the number of exon blocks for a single bait mapped to a single locus. `contig_4562` with 35 lines simply spans 35 exons, which makes sense for your longest baits in the 4,000–6,000 bp range. A 5,000 bp CDS from a large gene spanning many short exons will produce exactly this pattern. This is correct and expected behaviour from `-split`.*

*Two small things to note though. First, your bait count is 1,218 in the BED but you have 1,250 loci in the bait FASTA. That means **32 baits produced no primary alignment** and were silently dropped. Worth identifying them with the script *`script_finding_missing_baits.sh`*.*



## 3. Checking what are the baits that are not mapping with the B. masoniana genome

Eventually, the result is:

wc -l baits_unmapped.txt

21 baits_unmapped.txt, but checking on the baits that are mapped in the file baits_mapped.txt, we find 1,218 baits correctly mapped to B.masoniana, which is the majority of conserved loci between B. luzhaiensis and B. masoniana. We can now proceed to include the annotations in the file.
