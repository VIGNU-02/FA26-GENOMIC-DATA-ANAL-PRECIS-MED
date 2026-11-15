# Assignment 2 | miRNA Differential Expression in SARS-CoV-2 Infected Respiratory Cells

**Vignesh Gutta** | INFO-B636 Genomic Data Analysis for Precision Medicine | 11/15/2026

Twelve small-RNA-seq libraries from cultivated human respiratory cells (BioProject PRJNA901149) were processed to find microRNAs that change on SARS-CoV-2 infection, and to find what the affected miRNAs regulate.

**Human genome build: GRCh38, Ensembl release 115. miRNA annotation: miRBase v23.1.**

## Design

| | 24 h | 72 h |
|---|---|---|
| Mock | SRR22269881, 882, 883 | SRR22269875, 876, 877 |
| SARS-CoV-2 | SRR22269878, 879, 880 | SRR22269872, 873, 874 |

Two comparisons: Mock vs SARS-CoV-2 across both time points (`~ time + condition`), and SARS-CoV-2 at 24 h vs 72 h (`~ time`, infected samples only).

## Pipeline

```
SRA download  ->  FastQC  ->  Trim Galore  ->  HISAT2 (GRCh38)  ->  featureCounts (miRBase)
                                                                             |
                              GO enrichment  <-  DESeq2  <-------------------+
```

All alignment and counting ran on IU Big Red 200 under Slurm (account c02374); the R steps ran in the same conda environment.

## Results

| | Mock vs SARS | SARS 24 h vs 72 h |
|---|---|---|
| miRNAs tested | 527 | 527 |
| Significant (padj < 0.05) | 28 (11 up, 17 down) | 49 (22 up, 27 down) |
| Strongest | miR-1246-5p, +1.80, padj 1.6e-10 | miR-200c-3p, +1.21, padj 1.7e-05 |

Five let-7 family members (let-7a-5p, let-7f-5p, let-7g-5p, let-7i-5p, miR-98-5p) all go down on infection, and the miR-26, miR-199, miR-374, miR-550 and miR-941 families each move coherently as well.

GO enrichment on the validated targets of the up-regulated miRNAs, run against a miRNA-target background rather than the whole genome, returns 13 terms, dominated by autophagy and nutrient sensing.

Full numbers, interpretation and the limitations of this dataset are in **[REPORT.md](REPORT.md)**.

## Reproducing this

Every command, software version, download URL and script, in run order, is in **[README.txt](README.txt)**.

Two things in that file are worth flagging up front, because both will silently produce wrong results if skipped:

- The adapter is the Illumina TruSeq small-RNA adapter `AGATCGGAAGAGC` and has to be passed explicitly. Trim Galore's `--small_rna` preset assumes a different one and removed adapter from only 22.5% of reads here.
- miRBase names chromosomes `chr1`, the Ensembl genome names them `1`. The prefix has to be stripped (`sed 's/^chr//'`) or featureCounts matches nothing and every count comes out zero.

## Layout

Folders follow the pipeline in order. Each one holds the script for that step and the output it produced.

```
README.txt                          all commands, versions, GRCh38 build, run order
README.md                           this file
REPORT.md                           results and interpretation

01_trimming/
    trim_array.sh                   Trim Galore, Slurm array over 12 samples

02_alignment/
    build_index.sh                  hisat2-build on GRCh38
    align_array.sh                  HISAT2 + samtools sort, Slurm array
    alignment_rates.txt             overall alignment rate per sample

03_counting/
    count.sh                        featureCounts against miRBase v23.1
    count_ensembl.sh                featureCounts against the Ensembl GTF (comparison only)
    metadata.csv                    sample sheet: condition and time point
    mirna_counts_clean.tsv          raw count matrix, 2,783 miRNAs x 12 samples
    mirna_counts.txt.summary        featureCounts assignment summary

04_differential_expression/
    deseq2_analysis.R               DESeq2, both comparisons, PCA
    volcano_plots.R                 volcano plots
    DEG_Mock_vs_SARS_full.csv               all 527 tested
    DEG_Mock_vs_SARS_significant.csv        28 at padj < 0.05
    DEG_SARS_24h_vs_72h_full.csv            all 527 tested
    DEG_SARS_24h_vs_72h_significant.csv     49 at padj < 0.05
    miRNA_normalized_expression.csv         DESeq2 normalized counts
    PCA_plot.pdf
    volcano_Mock_vs_SARS.png
    volcano_SARS_24h_vs_72h.png

05_go_enrichment/
    go_final.R                      GO enrichment with corrected background
    GOfinal_Mock_vs_SARS_UP.csv             13 terms
    GOfinal_Mock_vs_SARS_DOWN.csv           749 terms
    GOfinal_SARS_24h_vs_72h_UP.csv          307 terms
    GOfinal_SARS_24h_vs_72h_DOWN.csv        573 terms
```

The SRA downloads, FastQC output, trimmed FASTQ files, BAM files and the GRCh38 index are not in the repository; they are too large. `README.txt` has the commands that regenerate all of them.
