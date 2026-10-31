# Midterm Report | Oxford Nanopore Direct RNA Sequencing Analysis

**Programmer:** Vignesh Gutta (vgutta@iu.edu)
**Date:** 10/31/2026
**System:** IU Big Red 200 (Slurm account c02374)
**Dataset:** `rna004_hsapiens_r10k.pod5` from Zenodo record 15853348. 10,000 human direct RNA reads, RNA004 chemistry, Universal Human Reference RNA (UHR), poly-A selected.

---

## 1. Basecalling

### Model and justification
Model: `rna004_130bps_sup@v6.0.0`, auto-selected by Dorado when I passed the `sup` alias.

I chose `sup` over `fast` and `hac` because:
- The dataset is only 10,000 reads, so runtime was not a constraint.
- On an NVIDIA A100-SXM4-40GB the sup model finished in 35 seconds.
- Higher basecalling accuracy gives more reads that align in Step 3 and cleaner signal-to-base assignment in Step 4.
- ONT recommends sup for production analysis unless throughput is the priority.

### Parameters
- `--emit-moves`: saves the move table, which Squigualiser needs in Step 4. Without it, Step 4 cannot run.
- `--device cuda:0`: run on the GPU.

### Command
```bash
dorado basecaller sup rna004_hsapiens_r10k.pod5 \
  --emit-moves --device cuda:0 > basecall/basecalls.bam
```

### Results
- 9,997 reads basecalled, 49 filtered by Dorado's internal QC.
- 53 reads were split into chimeric sub-reads, so the final BAM holds 10,050 records.
- Runtime: 35 seconds on one A100.
- Output: 18 MB BAM.

---

## 2. Quality Assessment

NanoPlot 1.48.0 run on the unaligned BAM:

| Metric | Value |
|---|---|
| Reads | 10,050 |
| Total bases | 11.2 Mb |
| Mean read length | 1,115 bp |
| Median read length | 828 bp |
| N50 | 1,526 bp |
| Longest read | 49,109 bp |
| Mean quality | Q14.5 |
| Median quality | Q18.7 |
| Reads >Q10 | 9,484 (94.4%) |
| Reads >Q15 | 8,506 (84.6%) |
| Reads >Q20 | 2,967 (29.5%) |

### Interpretation
Median quality is Q18.7, roughly one error per 75 bases, which sits at the high end for direct RNA and reflects the RNA004 chemistry combined with sup basecalling. 85% of reads pass Q15.

The length distribution fits full-length direct-RNA capture: median 828 bp, N50 1,526 bp. The longest reads (over 20 kb) carry quality scores of Q1 to Q2, so they are chimeras or stall artifacts rather than real molecules. The highest-quality reads (Q26 to Q29) are 275 to 1,423 bp long, the range where direct RNA performs best.

Total yield of 11.2 Mb is in line with a 10,000-read subset drawn from a larger run.

### Command
```bash
NanoPlot --ubam basecall/basecalls.bam -o nanoplot --tsv_stats --raw
```

I used `--ubam` rather than `--bam` because the BAM is unaligned at this stage. Report: `NanoPlot-report.html`.

---

## 3. Alignment

### Tool and justification
Aligner: **minimap2 2.31**, `-ax map-ont` preset.

minimap2 is the standard long-read aligner for Oxford Nanopore data. Its seed-chain-align strategy (minimizers, then chaining, then base-level alignment) handles noisy long reads efficiently, and ONT recommends it for direct RNA. I used it rather than a short-read aligner like BWA because short-read aligners assume low error rates and short fragments, neither of which holds for nanopore RNA.

### Reference
Ensembl *Homo sapiens* GRCh38 transcriptome, release 112. I combined `cdna.all.fa` (coding) and `ncrna.fa` (noncoding) into one file of 277,062 transcripts. The cDNA file alone does not contain mitochondrial rRNAs such as MT-RNR2, the Step 4 target, so the noncoding file is required.

### Parameters
- `-ax map-ont`: Nanopore preset, SAM output.
- `-uf`: forward strand only. Direct RNA is strand-specific because the molecule's 3′ end enters the pore first.
- `-k 14`: k-mer length ONT recommends for direct RNA.
- `-y`: carry input SAM tags into the output, which keeps Dorado's move table (`mv:`) for Step 4.
- `--MD`: emit the MD tag for per-base mismatches.
- `-t 8`: 8 threads.

### Command
```bash
samtools fastq -T "*" basecall/basecalls.bam > basecall/basecalls.fastq

minimap2 -ax map-ont -uf -k 14 -t 8 -y --MD \
  reference/Homo_sapiens.GRCh38.transcriptome.fa \
  basecall/basecalls.fastq | samtools sort -o align/aligned.bam -
samtools index align/aligned.bam
```

I converted BAM to FASTQ first because minimap2 reads FASTA/FASTQ, not BAM. `samtools fastq -T "*"` carries the move-table tags through the conversion.

### Results

| Metric | Value |
|---|---|
| Primary reads | 10,050 |
| Primary mapped | 9,562 (95.14%) |
| Secondary alignments | 24,254 |
| Supplementary alignments | 533 |
| Total alignment records | 34,837 |
| Reads on MT-RNR2 (ENST00000387347.2) | 183 |

### Top 10 transcripts by read count

| Count | Gene | Function |
|---|---|---|
| 183 | MT-RNR2 | 16S mitochondrial rRNA |
| 142 | MT-CO1 | Cytochrome c oxidase subunit 1 |
| 130 | MT-CO3 | Cytochrome c oxidase subunit 3 |
| 105 | MT-ND4 | NADH dehydrogenase subunit 4 |
| 101 | MT-TS1 | Mitochondrial tRNA (serine) |
| 90 | MT-ATP6 | ATP synthase subunit 6 |
| 88 | MT-CO2 | Cytochrome c oxidase subunit 2 |
| 74 | MT-RNR1 | 12S mitochondrial rRNA |
| 69 | MT-ATP8 | ATP synthase subunit 8 |
| 60 | MT-CYB | Cytochrome b |

### Interpretation
95.14% of primary reads mapped, which tells me the sample is clean human RNA and the basecalls are accurate enough to place.

The ~24,000 secondary alignments come from transcript isoforms: a read covering an exon shared across several isoforms of one gene gets reported once per isoform. Resolving which isoform produced each read would need isoform-level quantification (for example NanoCount or Salmon).

All of the top 10 transcripts are mitochondrial. A cell carries hundreds to thousands of mitochondria, and their long transcripts (the MT-RNR rRNAs, the MT-CO/ND/ATP mRNAs, the MT-tRNAs) are polyadenylated, so they survive the poly-A selection. These 10 transcripts hold about 11% of mapped reads; the other 89% spread across thousands of nuclear transcripts at low per-transcript coverage, which is what an unamplified bulk direct-RNA sample looks like.

The 4.86% that did not map are likely chimeras Dorado did not split, very low-quality reads, and transcripts not yet in Ensembl.

---

## 4. Raw Signal Visualization

### Tool
Squigualiser 0.7.0, Option 2 (basecaller move table).

### Approach
Dorado's move table records which stretch of signal produced which base. Squigualiser uses it in three steps:

1. **reform**: convert the move table into a PAF that links signal samples to each read's bases.
2. **realign**: combine that PAF with the minimap2 BAM to map signal samples to reference positions.
3. **plot**: draw raw current (y-axis, picoamps) against signal sample index (x-axis), with reference bases labeled along the top.

I converted the POD5 to BLOW5 with `blue-crab` first, because Squigualiser's slow5 library reads BLOW5/SLOW5, not POD5.

### Commands
```bash
blue-crab p2s rna004_hsapiens_r10k.pod5 -o rna004_hsapiens_r10k.blow5

squigualiser reform --sig_move_offset 0 --kmer_length 5 \
  -b basecall/basecalls.bam -c -o squig/reform.paf

squigualiser realign --bam align/aligned.bam \
  --paf squig/reform.paf -o squig/realign.bam
samtools index squig/realign.bam

squigualiser plot \
  --file reference/Homo_sapiens.GRCh38.transcriptome.fa \
  --slow5 rna004_hsapiens_r10k.blow5 \
  --alignment squig/realign.bam \
  --region ENST00000387347.2:1360-1400 \
  --output_dir squig/plot \
  --tag_name "MT-RNR2_1360-1400" \
  --plot_limit 10 --rna
```

### Parameters
- `--sig_move_offset 0`, `--kmer_length 5`: the RNA004 pore holds 5 bases at a time with no offset between move entries and signal. Squigualiser detected a stride of 6 signal samples per move entry.
- `--rna`: RNA data, read 3′ to 5′.
- `--plot_limit 10`: up to 10 plots, one per read covering the region.

### Results
183 reads mapped to MT-RNR2. Of those, 10 cover the full 1360-1400 window and were plotted as interactive Bokeh HTML files (`midterm_plots/*.html`).

### What the plot shows
The x-axis is the signal sample index, in the order the pore recorded measurements. The y-axis is raw ionic current in picoamps. The top strip shows the MT-RNR2 reference from position 1360 to 1400, one label per base. Under each base sits the current recorded while that base and its surrounding 5-mer were inside the pore.

Different bases give different current levels because the narrow part of the pore blocks ionic flow by an amount set by the specific 5-mer sitting in it. The strand reads 3′ to 5′ because direct-RNA sequencing starts at the poly-A tail.

### Interpretation
Across the 10 reads, the current level and shape under each labeled base line up. The pore produces a repeatable signature for the same base context across different molecules, and Dorado's basecalls track that signal. A position where every read showed an unexpected current shape would point to an RNA modification such as m6A, since modifications shift the current.

---

## 5. Precision Medicine Application

**Isoform-resolved transcriptomics for cancer therapy selection.**

Short-read RNA-seq produces 100 to 150 bp fragments that cover only part of a transcript. When a gene has several splice isoforms that share exons, a short read cannot say which full-length isoform it came from, so isoform-level counts rely on statistical inference and often collapse to the gene. A long direct-RNA read covers the whole mRNA in one pass and names its isoform directly.

Three cases where that matters:

1. **Fusion transcripts in cancer.** Oncogenic fusions such as BCR-ABL1 in CML or EML4-ALK in lung cancer create chimeric transcripts. A short read across the fusion junction is rare and often ambiguous; a long read spans the whole fusion and gives the exact breakpoint and reading frame, which decides whether a targeted drug like imatinib or crizotinib is likely to work.

2. **Isoform switching.** Splice variants of the androgen receptor (AR-V7) in prostate cancer predict resistance to abiraterone and enzalutamide. Long reads quantify AR-V7 against full-length AR without a variant-specific assay, so a novel variant is not missed.

3. **RNA modifications as biomarkers.** Direct RNA sequencing keeps chemical modifications (m6A, pseudouridine, inosine) in the signal, because the native molecule passes through the pore. m6A changes on specific transcripts have been linked to AML aggressiveness and could help stratify patients. Short-read RNA-seq converts RNA to cDNA and erases these marks.

In each case the long read carries information a short read cannot: the full isoform, the exact breakpoint, or the native modification.
