# summary

Use **phylogeny-aware** and **comparative** approaches: codon-models to detect selection (branch/site and site tests), conservation-based deleteriousness predictions, protein-level functional predictions and structural mapping, and convergence / phenotype-association tests across species — then combine lines of evidence to prioritise candidates.

---

# 1) Detect positive (adaptive) selection with codon/phylogenetic models

Because you have a phylogeny and target loci across many species, codon models are the standard choice.

Practical tests to run (per-gene / per-alignment):

* **Branch-level**: aBSREL — tests whether particular branches experienced episodic diversifying selection (good to find lineages with bursts of positive selection).
* **Gene / gene-wide**: BUSTED — test for evidence of positive selection anywhere in the gene. ([help.datamonkey.org][2])
  http://help.datamonkey.org/methods/busted.html
* **Site-level**: MEME — finds sites that experienced episodic positive selection on a subset of branches (more powerful than classic site tests when selection is episodic). Use FEL/FUBAR for pervasive selection as complements. ([stevenweaver.github.io][3])

Tooling: **HyPhy** (command-line or Datamonkey) implements aBSREL, BUSTED, MEME, FEL, FUBAR and is well documented. ([veg.github.io][1])

Notes:

* Use *codon* alignments (translate, align proteins, convert back to codons with PAL2NAL).
* Remove paralogs and poor alignments; trim misaligned ends (trimAl/Gblocks) to avoid false positives.
* Correct for multiple testing across many genes.

---

# 2) Detect relaxed vs intensified selection

If you want to detect genes under **relaxed** constraint (which may let deleterious alleles accumulate) or intensified selection, use **RELAX** (HyPhy) to quantify relaxation/intensification of selection regimes on specified sets of branches. This is useful when adaptation may involve relaxation followed by selection. ([veg.github.io][1])

---

# 3) Per-site evolutionary conservation → predict deleteriousness

Positions that are highly conserved across the phylogeny are more likely to be functionally important; non-synonymous changes at those sites are candidate **deleterious** variants.

Options:

* **Rate4Site / ConSurf** — produce per-site conservation scores from a protein MSA + tree (gives relative conservation rates). Use these scores to flag mutations at conserved sites. ([Tel Aviv University][4])
* **phylogenetic GERP-style**: if you have deep species sampling you can estimate constraint per site (but classical GERP requires dense alignments & neutral model — Rate4Site/ConSurf are simpler).

---

# 4) Functional effect predictors (protein level)

For coding changes you can run in-silico effect predictors that use sequence conservation / homology:

* **PROVEAN, SIFT** — produce scores predicting whether an amino-acid substitution is likely deleterious (PROVEAN works well for non-human species if you provide homologs). ([J. Craig Venter Institute][5])
* **PolyPhen-2** is human-centric (requires building a species database) — usually harder for non-model organisms but still possible.

Combine these with conservation scores above: change at a highly conserved site *and* a damaging PROVEAN/SIFT call = stronger evidence of deleterious effect.

---

# 5) Structural mapping / AlphaFold

Map candidate amino-acid changes onto protein structures (experimental or predicted by AlphaFold). Structural context (active site, binding interface, buried core) is strong evidence that a substitution is functionally important or deleterious.

---

# 6) Detect convergent / phenotype-associated substitutions

If you have phenotype or ecological labels for species (e.g. habitat, trait state), look for **convergent** amino-acid changes associated with the trait:

* **ConDor** — two-component approach (emergence + correlation) to find convergent substitutions in protein alignments. Good for large protein alignments. ([condor.pasteur.cloud][6])
* **PCOC / phyloConverge** — methods for detecting convergent shifts in amino-acid patterns or rates correlated with phenotype. Use BayesTraits or phylogenetic logistic/continuous models for trait–substitution correlation as well. ([Royal Society Publishing][7])

Caveat: convergent signal can be confounded by mutational biases or shared ancestry — always verify with simulations or null models.

---

# 7) Use ancestral reconstruction and branch patterns

Reconstruct ancestral sequences (e.g. using PAML or HyPhy) to:

* infer *where* substitutions happened (terminal vs internal branches),
* determine whether a change is unique to a lineage (possible deleterious private substitution) or fixed across a clade (candidate adaptation).

---

# 8) Integrative prioritisation (best practice)

Because each method has false positives, combine evidence to prioritise candidates:

High-priority *adaptive* candidate example:

* gene shows significant aBSREL/BUSTED signal on branches of interest **AND**
* specific site(s) identified by MEME **AND**
* site is not hypervariable (conserved by Rate4Site) **AND**
* substitutions correlate with phenotype/convergence test.

High-priority *deleterious* candidate example:

* nonsynonymous change occurs at a very conserved site **AND**
* predicted damaging by PROVEAN/SIFT **AND**
* occurs on terminal branch (private) or shows signatures of relaxed selection at the gene (RELAX).

Then do GO enrichment and pathway checks to see if candidates cluster in functions (helps interpret).

---

# 9) Practical pipeline / commands (very brief)

1. Orthology & alignment: identify ortholog loci → translate → protein MSA (MAFFT) → back-translate (PAL2NAL) → trim (trimAl).
2. Tree per gene or use your concatenated tree (make sure gene trees don’t conflict badly).
3. HyPhy analyses:

   * MEME: `hyphy MEME --alignment gene_codon.fasta --tree gene_tree.nwk`
   * aBSREL: `hyphy absrel --alignment gene_codon.fasta --tree gene_tree.nwk`
   * RELAX: `hyphy RELAX --alignment ... --tree ... --test`
     (See HyPhy docs for options). ([stevenweaver.github.io][8])
4. Conservation: `rate4site` or use ConSurf webserver. ([Tel Aviv University][4])
5. PROVEAN: run against protein sequences with homologs or local database. ([J. Craig Venter Institute][5])
6. Convergence: run ConDor/PCOC and/or phyloConverge if you have phenotypes. ([condor.pasteur.cloud][6])

(If you want, I can paste a short example script for steps 1–3 for one gene.)

---

# 10) Important caveats & data QC

* **Target capture**: loci are not random — enrichment biases and missing data can affect tests (especially if alignments are shallow for some loci). Interpret carefully.
* **Sample sizes per species**: while tests above use substitution patterns across the tree, very shallow sampling per species still reduces power for some inferences (e.g. polymorphism vs fixed).
* **Paralogs / assembly errors**: these create false positives — filter stringently for orthology.
* **Recombination** within loci can bias codon models — screen with GARD (HyPhy) and partition if needed. ([hyphy.org][9])

---

# Final notes & offer

If you’d like I can:

* (A) give you a ready-to-run shell/Python pipeline for a single locus that does the align→PAL2NAL→trim→HyPhy MEME/aBSREL→Rate4Site→PROVEAN steps; **or**
* (B) look at one or two example loci/alignments (you uploaded the tree at `/mnt/data/concat_finaltree.treefile`) and produce a worked example (I’ll need the codon alignments or raw FASTAs per locus + a mapping of individuals→species/phenotypes).

Tell me which you prefer and I’ll produce the exact scripts/commands and example outputs — or I can just paste a minimal pipeline here now.

[1]: https://veg.github.io/hyphy-site/tutorials/CL-prompt-tutorial/?utm_source=chatgpt.com
[2]: https://help.datamonkey.org/methods/meme?utm_source=chatgpt.com
[3]: https://stevenweaver.github.io/hyphy-site/methods/selection-methods/?utm_source=chatgpt.com
[4]: https://www.tau.ac.il/~itaymay/cp/rate4site.html?utm_source=chatgpt.com
[5]: https://www.jcvi.org/research/provean?utm_source=chatgpt.com
[6]: https://condor.pasteur.cloud/?utm_source=chatgpt.com
[7]: https://royalsocietypublishing.org/doi/10.1098/rstb.2018.0234?utm_source=chatgpt.com
[8]: https://stevenweaver.github.io/hyphy-site/getting-started/?utm_source=chatgpt.com
[9]: https://hyphy.org/resources/tutorial-2017.pdf?utm_source=chatgpt.com
