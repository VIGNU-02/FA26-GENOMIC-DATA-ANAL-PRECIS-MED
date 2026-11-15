# Assignment 2 Report | miRNA Expression in SARS-CoV-2 Infected Respiratory Cells

**Programmer:** Vignesh Gutta
**Course:** INFO-B636, Genomic Data Analysis for Precision Medicine
**Date:** 11/15/2026
**System:** IU Big Red 200 (Slurm account c02374)
**Data:** BioProject PRJNA901149, 12 small-RNA-seq libraries from cultivated human respiratory cells
**Genome build:** GRCh38 (Ensembl release 115)
**miRNA annotation:** miRBase v23.1 (`hsa.gff3`, GRCh38 coordinates)

---

## 1. Experimental design

Twelve samples in a 2 x 2 design, three biological replicates per group:

| | 24 h | 72 h |
|---|---|---|
| Mock (uninfected) | SRR22269881, 882, 883 | SRR22269875, 876, 877 |
| SARS-CoV-2 infected | SRR22269878, 879, 880 | SRR22269872, 873, 874 |

Two comparisons were requested:

1. **Mock vs SARS-CoV-2**, irrespective of time point. Run as a DESeq2 design of `~ time + condition`, which controls for the time effect and isolates infection.
2. **SARS-CoV-2 at 24 h vs 72 h**, using the six infected samples only, design `~ time`.

---

## 2. Data processing summary

The libraries are small-RNA-seq. All raw reads were exactly 50 bp, with the Illumina TruSeq adapter (`AGATCGGAAGAGC`) beginning partway through most reads, because mature miRNAs are only about 22 nt.

A first trimming run using Trim Galore's `--small_rna` preset removed adapter from only 22.5% of reads, because that preset assumes a different adapter (`TGGAATTCTCGG`). Specifying the adapter explicitly raised adapter detection to 94.7% and dropped residual adapter in the trimmed output from 3,879,925 occurrences to zero. After trimming, the modal read length was 22 nt, which matches the expected mature miRNA size.

Alignment rates to GRCh38 ran from 23.7% to 72.9%, median around 61%. Rates in this range are expected for small-RNA libraries, which carry a lot of non-miRNA material alongside the miRNAs: tRNA and rRNA fragments, degradation products. SRR22269883 was the outlier at 23.7%.

Counting against miRBase assigned 368,752 to 4,330,206 reads per sample to mature miRNA features.

| Metric | Value |
|---|---|
| miRNAs in miRBase catalogue | 2,783 |
| Detected (>0 reads in any sample) | 1,645 |
| Expressed at >= 100 total reads | 485 |
| Retained for DE testing (>= 10 counts in >= 3 samples) | 527 |

---

## 3. Choice of annotation: miRBase vs Ensembl

The assignment asks for expression of microRNA genes, so reads were counted against the miRBase mature-miRNA annotation rather than the full Ensembl gene GTF. Both were run on the same alignments to document the difference:

| | miRBase `hsa.gff3` | Ensembl GRCh38.115 GTF |
|---|---|---|
| Features counted | 2,783 (all mature miRNAs) | 78,899 genes (all biotypes) |
| Reads assigned, SRR22269872 | 1,816,318 | 5,734,791 |
| Of which are miRNA genes | 100% | 1,879 / 78,899 = 2.4% |

The full Ensembl GTF assigns roughly three times more reads, but 97.6% of the features it measures are not miRNAs. In a size-selected small-RNA library those extra assignments are mostly fragments of lncRNAs, snoRNAs, tRNAs and protein-coding transcripts, not real expression of those genes. The miRBase annotation was used for everything downstream.

---

## 4. Differential expression

DESeq2 (median-of-ratios normalization, Wald test, Benjamini-Hochberg correction) on the 527 retained miRNAs.

### Threshold choice

Two thresholds were examined for comparison 1:

| Criterion | DEGs |
|---|---|
| padj < 0.05 | 28 |
| padj < 0.05 and abs(log2FC) > 1 | 4 |

The common abs(log2FC) > 1 (2-fold) cutoff comes from mRNA studies. miRNAs tune their targets rather than switch them, so their changes are smaller, and applying that cutoff here threw out 24 miRNAs that were statistically solid, among them miR-877-5p at log2FC = 0.71 and padj = 6.1e-05. padj < 0.05 was used as the primary criterion, with fold change reported alongside so the size of each change stays visible.

### Comparison 1: Mock vs SARS-CoV-2

28 differentially expressed miRNAs: 11 up, 17 down in infection.

Most significant:

| miRNA | log2FC | baseMean | padj |
|---|---|---|---|
| hsa-miR-1246-5p | +1.803 | 95 | 1.64e-10 |
| hsa-miR-877-5p | +0.708 | 733 | 6.15e-05 |
| hsa-miR-4728-5p | +1.298 | 52 | 6.93e-05 |
| hsa-miR-4521-5p | -0.830 | 93 | 3.04e-03 |
| hsa-miR-21-3p | -1.340 | 85 | 5.47e-03 |
| hsa-let-7g-5p | -0.907 | 53,828 | 5.47e-03 |
| hsa-miR-574-3p | +0.704 | 516 | 5.47e-03 |
| hsa-miR-98-5p | -0.841 | 1,398 | 5.47e-03 |

Full list: `04_differential_expression/DEG_Mock_vs_SARS_significant.csv`. Volcano plot: `04_differential_expression/volcano_Mock_vs_SARS.png`.

Related miRNAs move together:

| Family | Members affected | Direction |
|---|---|---|
| let-7 | let-7a-5p, let-7f-5p, let-7g-5p, let-7i-5p, miR-98-5p | all down |
| miR-26 | miR-26a-5p, miR-26b-5p | both down |
| miR-199 | miR-199a-3p, miR-199b-3p | both down |
| miR-374 | miR-374b-5p, miR-374c-3p | both down |
| miR-550 | miR-550a-5p, miR-550b-3p | both up |
| miR-941 | miR-941-5p, miR-941-3p | both up |

Five let-7 family members go down on infection, and these are among the most abundant miRNAs in the dataset (baseMean 53,828 to 143,651), so their expression is measured with high precision. Random variation does not move whole sequence-related families in one direction, which is the main internal check this result passes.

let-7 family members repress pro-inflammatory transcripts including IL-6, so coordinated let-7 loss de-represses inflammatory signaling, a pathway of considerable interest in COVID-19.

miR-1246, the most significant miRNA here, has been reported in the literature as elevated in SARS-CoV-2 infection, which supports the pipeline from outside the dataset.

### Comparison 2: SARS-CoV-2 24 h vs 72 h

49 differentially expressed miRNAs: 22 up, 27 down at 72 h relative to 24 h.

Most significant:

| miRNA | log2FC (72h vs 24h) | baseMean | padj |
|---|---|---|---|
| hsa-miR-200c-3p | +1.209 | 20,166 | 1.66e-05 |
| hsa-miR-4521-5p | -1.664 | 69 | 2.73e-04 |
| hsa-miR-224-5p | +1.016 | 3,686 | 6.12e-04 |
| hsa-miR-148b-3p | -1.453 | 4,276 | 1.01e-03 |
| hsa-miR-148a-3p | -1.508 | 20,215 | 1.05e-03 |
| hsa-miR-99b-5p | +1.102 | 16,650 | 1.60e-03 |
| hsa-miR-199a-3p | -1.150 | 3,405 | 1.93e-03 |

Full list: `04_differential_expression/DEG_SARS_24h_vs_72h_significant.csv`. Volcano plot: `04_differential_expression/volcano_SARS_24h_vs_72h.png`.

Family coherence shows up again: miR-148a-3p and miR-148b-3p both fall, as do miR-199a-3p and miR-199b-3p, the pair that also changed in comparison 1.

More miRNAs change with time inside infection (49) than between infected and uninfected cells overall (28), and the changes are larger. The host miRNA response develops over the course of infection, with the 72 h state well separated from 24 h, rather than settling into one state as soon as the cells are infected.

---

## 5. GO enrichment

### Method, and why the first attempts were discarded

miRNAs regulate rather than act, and they carry sparse direct GO annotation, so enrichment was run on their target genes. Experimentally validated miRNA-target interactions came from miRTarBase through the multiMiR package (database v2.4.0), restricted to strongly supported interactions ("Functional MTI", weak evidence excluded). GO Biological Process enrichment was then run with clusterProfiler.

Three problems turned up and were corrected.

**Direction was ignored at first.** An increase in a miRNA suppresses its targets while a decrease de-represses them, so analyzing all DEGs together mixes opposing effects. Up- and down-regulated miRNAs were separated.

**The target set was too large.** All validated interactions for all 28 DEG miRNAs gave 5,349 target genes, about a quarter of the protein-coding genome, and 1,963 "enriched" terms. A query set that large carries enough statistical power that nearly every category clears significance, so nothing in the output can be read as specific. Restricting to strongly validated interactions cut this down a lot.

**The background was wrong.** `enrichGO` defaults to the whole genome as background. miRNA targets are not a random sample of genes; they skew toward well-studied, heavily annotated genes such as cancer and developmental regulators. That bias was visible in the output, with "epithelial cell proliferation" and "gland development" turning up as top terms for three different input sets. The background was redefined as the union of strongly validated targets of all 527 expressed miRNAs (2,454 genes). The question then becomes whether, among genes that miRNAs regulate in this system, the targets of the differentially expressed miRNAs stand out.

Effect of the background correction:

| Gene set | Terms, genome background | Terms, miRNA-target background |
|---|---|---|
| Mock vs SARS, up | 1,963 | 13 |
| Mock vs SARS, down | 2,003 | 749 |
| 24 h vs 72 h, up | 1,674 | 307 |
| 24 h vs 72 h, down | 1,741 | 573 |

### Result: Mock vs SARS, up-regulated miRNAs

The 11 up-regulated miRNAs map to 34 strongly validated targets and give 13 enriched GO BP terms:

| GO term | Genes | padj |
|---|---|---|
| autophagy | 10/33 | 2.85e-02 |
| process utilizing autophagic mechanism | 10/33 | 2.85e-02 |
| regulation of autophagy | 9/33 | 2.85e-02 |
| response to lipid | 14/33 | 3.33e-02 |
| head development | 12/33 | 2.85e-02 |
| response to nutrient levels | 10/33 | 3.33e-02 |
| cellular response to nutrient levels | 8/33 | 2.85e-02 |
| fat cell differentiation | 7/33 | 3.33e-02 |
| positive regulation of miRNA transcription | 6/33 | 2.85e-02 |
| regulation of miRNA transcription | 6/33 | 2.85e-02 |
| miRNA transcription | 6/33 | 2.85e-02 |
| positive regulation of miRNA metabolic process | 6/33 | 2.85e-02 |
| regulation of miRNA metabolic process | 6/33 | 3.33e-02 |

Three of the thirteen terms are autophagy-related, and nutrient-sensing terms, coupled to autophagy through mTOR signaling, account for two more. SARS-CoV-2 interacts extensively with host autophagy: it uses double-membrane vesicles derived from the autophagic machinery for replication and interferes with autophagosome-lysosome fusion. miRNAs that rise on infection targeting autophagy regulators fits that biology.

"miRNA transcription" and "miRNA metabolic process" terms also recur, so part of the response runs through the miRNA machinery regulating itself.

### Result: the other three sets

| Set | Terms | Leading themes |
|---|---|---|
| Mock vs SARS, down | 749 | negative regulation of gene expression; positive regulation of cell proliferation; negative regulation of programmed cell death; PI3K/AKT signaling |
| 24 h vs 72 h, up | 307 | negative regulation of programmed cell death; negative regulation of apoptotic process; circulatory system development; hemopoiesis |
| 24 h vs 72 h, down | 573 | circulatory system development; tube development; regulation of cell migration; cell proliferation |

Apoptosis regulation leads the 24 h vs 72 h up-regulated set, which fits a cell-fate decision shifting as infection runs on. These three sets each keep several hundred terms, so individual terms are indicative rather than conclusive. Top-ranked terms are given here and the full tables are in the repository.

All four GO result tables: `05_go_enrichment/GOfinal_*.csv`.

---

## 6. Quality assessment and limitations

**Sequencing depth is uneven.** Reads assigned to miRNAs ran from 368,752 (SRR22269883) to 4,330,206 (SRR22269876), a 12-fold spread, with DESeq2 size factors from 0.234 to 2.187. DESeq2 normalizes for this, but scaling a shallow library up amplifies its noise.

**Depth is partly confounded with time point.** Mean library size was about 1.42 M at 24 h and 2.46 M at 72 h, a 1.7-fold difference, and that difference falls across the groups compared in comparison 2. Depth was not meaningfully confounded with infection status (Mock about 1.77 M, SARS about 2.11 M).

**PCA separates the groups only partly.** PC1 and PC2 explain 29% and 19% of variance, leaving over half unexplained. Mock 72 h and SARS 24 h separate cleanly on PC2, but the other groups overlap and within-group replicate spread is wide. SRR22269883, the 23.7%-alignment sample, sits away from its two Mock 24 h replicates. See `04_differential_expression/PCA_plot.pdf`.

**Effect sizes are modest.** Of the 28 significant miRNAs in comparison 1, only 4 exceed 2-fold change.

**Some methodological choices are open to challenge.** HISAT2 was used for alignment; bowtie is more conventional for very short small-RNA reads and might give higher mapping rates. Multi-mapping and multi-overlapping reads were counted (`-M -O`), which is necessary because miRNA families share sequence but inflates counts for those families. Targets came from miRTarBase alone; TarBase was dropped after it failed to respond, though it also contributes lower-confidence high-throughput data that would have enlarged the target set further.

**What holds up.** The coordinated behavior of the miRNA families, the very low p-value and literature support for miR-1246, and the specificity of the autophagy enrichment all survive the noise above. Results sitting near the padj = 0.05 boundary should be treated as provisional.

---

## 7. Conclusions

1. 28 miRNAs are differentially expressed between mock and SARS-CoV-2 infected respiratory cells (11 up, 17 down; padj < 0.05), the strongest being miR-1246-5p at +1.80 and padj = 1.6e-10.
2. Five members of the let-7 family go down together on infection, alongside coherent shifts in the miR-26, miR-199, miR-374, miR-550 and miR-941 families.
3. 49 miRNAs differ between 24 h and 72 h post-infection (22 up, 27 down), more changes and larger ones than between infected and uninfected cells overall.
4. GO enrichment on strongly validated targets of the up-regulated miRNAs, against a miRNA-target background, gives 13 terms dominated by autophagy and nutrient sensing, consistent with documented SARS-CoV-2 interactions with host autophagy.
5. The dataset carries real technical variability, listed in section 6. The conclusions above rest on the findings that are internally coherent and externally supported, not on borderline statistics.
