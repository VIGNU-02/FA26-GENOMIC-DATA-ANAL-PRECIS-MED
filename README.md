# Midterm | Oxford Nanopore Direct RNA Sequencing Analysis

Analysis of 10,000 human direct RNA reads (RNA004 chemistry) from Zenodo record 15853348. Covers basecalling, quality assessment, alignment, and raw signal visualization.

See `REPORT.md` for the written report and interpretation. See `README.txt` for reproduction steps and all commands.

## Repository contents

| File | Purpose |
|---|---|
| `REPORT.md` | Full written report (answers to all assignment questions) |
| `README.md` | This file |
| `README.txt` | Reproducibility: every command, version, parameter, reference used |
| `NanoPlot-report.html` | Quality assessment output (Step 2 deliverable) |
| `NanoStats.txt` | NanoPlot summary table |
| `midterm_plots/*.html` | Squigualiser signal plots for MT-RNR2 positions 1360-1400 (Step 4 deliverable; 10 reads) |
| `top10_transcripts.txt` | Most-represented transcripts after alignment |

## Headline results

- **Basecalling** (Dorado 2.1.2, `rna004_130bps_sup@v6.0.0`): 10,050 reads from 10,000 input (chimera splits accounted for), 35 seconds on an NVIDIA A100
- **Quality** (NanoPlot 1.48.0): median Q18.7, N50 1,526 bp, 85% of reads above Q15
- **Alignment** (minimap2 2.31, Ensembl GRCh38 release 112 cDNA + ncRNA, 277,062 transcripts): **95.14% primary mapped**, top 10 transcripts all mitochondrial
- **Signal plot** (Squigualiser 0.7.0, Option 2): 10 reads plotted across MT-RNR2 (ENST00000387347.2) positions 1360-1400
