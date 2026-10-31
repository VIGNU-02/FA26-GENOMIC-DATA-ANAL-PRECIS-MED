# Midterm Report | Oxford Nanopore Direct RNA Sequencing Analysis

**Programmer:** Vignesh Gutta (vgutta@iu.edu)
**Date:** 10/31/2026
**System:** IU Big Red 200 (Slurm account c02374)
**Dataset:** `rna004_hsapiens_r10k.pod5` from Zenodo record 15853348 — 10,000 human direct RNA reads, RNA004 chemistry, Universal Human Reference RNA (UHR), poly-A selected

---

## 1. Basecalling

### Model and justification
Model: `rna004_130bps_sup@v6.0.0` (auto-selected by Dorado when passing the `sup` alias)

I chose `sup` over `fast`/`hac` because:
- The dataset is small (10,000 reads), so runtime was never going to be a constraint
- Running on an NVIDIA A100-SXM4-40GB GPU, the sup model finishes in under a minute
- Higher accuracy improves downstream alignment success and signal-to-base matching for Step 4
- Current ONT guidance recommends sup for all production-quality analysis unless throughput is critical

### Parameters
- `--emit-moves`: required to save the move table, needed by Squigualiser (Step 4)
- `--device cuda:0`: use the GPU
- Default beam search and chunking

### Command
```bash
dorado basecaller sup rna004_hsapiens_r10k.pod5 \
  --emit-moves --device cuda:0 > basecall/basecalls.bam
```

### Results
- 9,997 reads basecalled, 49 reads filtered by Dorado's internal QC
- 53 reads split into chimeric sub-reads → final BAM contains 10,050 records
- Total runtime: 35 seconds on one A100
- Output: 18 MB BAM file

---

## 2. Quality Assessment

NanoPlot 1.48.0 was run on the unaligned BAM. Key metrics:

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
- **Quality is excellent.** Median Q18.7 corresponds to roughly 1 error per 75 bases; this is at the high end for direct RNA, reflecting the combined improvements of RNA004 chemistry and sup basecalling.
- **Length distribution is biologically reasonable.** Median 828 bp with N50 1,526 bp matches expected direct-RNA full-length transcript capture. The longest reads (>20 kb) have low quality (Q1–Q2) and are likely chimeras or stall artifacts, not real molecules.
- **High-quality reads sit at moderate length** (Q26-29 reads are 275–1,423 bp), which is the sweet spot for direct RNA.
- **Yield (11.2 Mb)** is appropriate for a 10,000-read subset drawn from a larger run.

### Command
```bash
NanoPlot --ubam basecall/basecalls.bam -o nanoplot --tsv_stats --raw
```

Report: `NanoPlot-report.html` (included in this repo).

---

## 3. Alignment

### Tool and justification
Aligner: **minimap2 2.31** with the `-ax map-ont` preset.

minimap2 is the standard long-read aligner for Oxford Nanopore data. It uses a seed-chain-align strategy (minimizers → chaining → base-level alignment) that handles noisy long reads efficiently. It is the aligner recommended by ONT in all official documentation and used in essentially every direct-RNA publication.

### Reference
**Ensembl Homo sapiens GRCh38 transcriptome, release 112.** I combined `cdna.all.fa` and `ncrna.fa` to form a complete transcriptome (277,062 transcripts). The coding-only cDNA file alone does not contain mitochondrial rRNAs such as MT-RNR2 (the target for Step 4), so including the noncoding file is essential.

### Parameters
- `-ax map-ont`: Nanopore preset; outputs SAM
- `-uf`: forward strand only (direct RNA is strand-specific; the molecule's 3′-end enters the pore first and the strand orientation is known)
- `-k 14`: k-mer length recommended by ONT for direct RNA
- `-y`: carry over all SAM tags from the input BAM into the aligned output — this preserves Dorado's move table (`mv:` tag), required for Step 4
- `--MD`: emit MD tag for per-base mismatch information
- `-t 8`: 8 CPU threads

### Command
```bash
samtools fastq -T "*" basecall/basecalls.bam > basecall/basecalls.fastq

minimap2 -ax map-ont -uf -k 14 -t 8 -y --MD \
  reference/Homo_sapiens.GRCh38.transcriptome.fa \
  basecall/basecalls.fastq | samtools sort -o align/aligned.bam -
samtools index align/aligned.bam
```

Note: I converted BAM → FASTQ first because minimap2 does not read BAM input natively (it expects FASTA/FASTQ). `samtools fastq -T "*"` carries the move table tags through the conversion.

### Results

| Metric | Value |
|---|---|
| Primary reads | 10,050 |
| **Primary mapped** | **9,562 (95.14%)** |
| Secondary alignments | 24,254 |
| Supplementary alignments | 533 |
| Total alignment records | 34,837 |
| Reads aligned to MT-RNR2 (ENST00000387347.2) | 183 |

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
- **95.14% primary alignment** is excellent for direct RNA. The high rate confirms that both sample quality and basecalling accuracy are good.
- **~24,000 secondary alignments** reflect transcript isoforms: a read sharing an exon across multiple transcripts of the same gene gets reported once per matching isoform. Downstream isoform-level quantification would be needed to resolve these.
- **Mitochondrial transcripts dominate the top hits** because mitochondria are present in hundreds-to-thousands of copies per cell and their long, polyadenylated transcripts (MT-RNR1/2 rRNAs, MT-CO/ND/ATP mRNAs, MT-tRNAs) survive the poly-A selection step. These top 10 account for only ~11% of all mapped reads; the remaining 89% are distributed across thousands of nuclear transcripts with low per-transcript coverage, as expected for an unamplified direct-RNA bulk sample.
- **The 4.86% unmapped reads** likely include chimeras that Dorado didn't split, very low-quality reads, and transcripts not yet annotated in Ensembl.

---

## 4. Raw Signal Visualization

### Tool
**Squigualiser 0.7.0**, using Option 2 (basecaller move table) from the documentation.

### Approach
Dorado's move table (saved via `--emit-moves` in Step 1) records which stretch of signal produced which base. Squigualiser uses this table in three sub-steps:

1. **reform**: translate Dorado's move table into a PAF file linking signal samples to each read's bases
2. **realign**: combine reform PAF with the minimap2 BAM to produce signal↔reference position mappings
3. **plot**: draw interactive HTML plots of the raw current (y-axis, picoamps) over signal samples (x-axis), with reference bases labeled on top

The POD5 was converted to BLOW5 using `blue-crab` because Squigualiser's slow5 library requires BLOW5/SLOW5 input, not POD5.

### Commands
```bash
# POD5 → BLOW5
blue-crab p2s rna004_hsapiens_r10k.pod5 -o rna004_hsapiens_r10k.blow5

# reform
squigualiser reform --sig_move_offset 0 --kmer_length 5 \
  -b basecall/basecalls.bam -c -o squig/reform.paf

# realign
squigualiser realign --bam align/aligned.bam \
  --paf squig/reform.paf -o squig/realign.bam
samtools index squig/realign.bam

# plot: ENST00000387347.2 positions 1360-1400
squigualiser plot \
  --file reference/Homo_sapiens.GRCh38.transcriptome.fa \
  --slow5 rna004_hsapiens_r10k.blow5 \
  --alignment squig/realign.bam \
  --region ENST00000387347.2:1360-1400 \
  --output_dir squig/plot \
  --tag_name "MT-RNR2_1360-1400" \
  --plot_limit 10 --rna
```

### Parameters explained
- `--sig_move_offset 0`, `--kmer_length 5`: standard RNA004 settings (the pore holds 5 bases at a time, no offset between move entries and signal samples). Squigualiser auto-detected a stride of 6 (6 signal samples per move entry).
- `--rna`: tells Squigualiser the data is RNA (direction, base alphabet)
- `--plot_limit 10`: produce up to 10 plots (one per read covering the region)

### Results
Of the 183 reads mapped to MT-RNR2, 10 reads cover the full 1360-1400 window and were plotted. Each plot is an interactive Bokeh HTML (`midterm_plots/*.html`).

### What the plot shows
- **X-axis:** signal sample index (chronological order of pore measurements)
- **Y-axis:** raw ionic current in picoamps (pA)
- **Top strip:** the reference sequence of MT-RNR2 from position 1360 to 1400, with each base labeled
- **Signal under each base:** the current trace recorded while that base (and its surrounding 5-mer context) was inside the pore
- Different bases produce distinct current levels because the physical constriction of the pore blocks ionic flow by amounts that depend on the specific 5-mer sitting inside it
- The strand is read 3′ → 5′ because direct-RNA sequencing starts at the poly-A tail

### Interpretation
Comparing multiple reads in the same region shows that signal patterns are consistent across molecules — the current levels under each labeled base match in shape and amplitude across different reads. This confirms:
1. Dorado's basecalls are a faithful interpretation of the raw signal
2. The pore produces reproducible current signatures for the same k-mer context across different molecules
3. Any position where all reads show an unexpected signal shape would be a candidate for an RNA modification (such as m6A), since modifications shift the current measurably

---

## 5. Precision Medicine Application

**Isoform-resolved transcriptomics in cancer therapy selection.**

Short-read RNA-seq (Illumina) produces 100-150 bp fragments that cover only part of a transcript. When a gene has multiple splice isoforms sharing many exons, short reads cannot determine which full-length isoform each fragment came from — isoform quantification relies on statistical inference and often collapses to the gene level.

Long-read direct RNA sequencing reads **entire mRNA molecules in one pass**, so each read definitively identifies its isoform. Three applications in precision medicine:

1. **Fusion transcripts in cancer.** Oncogenic fusions (e.g., BCR-ABL1 in CML, EML4-ALK in lung cancer) create novel chimeric transcripts. A short read spanning the fusion junction is often rare and ambiguous; a long read spans the whole fusion transcript and reveals the exact breakpoint and reading frame, which determines whether a targeted therapy (imatinib, crizotinib, etc.) is likely to work.

2. **Clinically relevant isoform switching.** Different cancer subtypes express different isoforms of the same gene. For example, splice variants of the androgen receptor (AR-V7) in prostate cancer predict resistance to abiraterone and enzalutamide. Long reads unambiguously quantify AR-V7 vs full-length AR; short reads require specialized assays that miss novel variants.

3. **RNA modifications as biomarkers.** Direct RNA sequencing uniquely preserves chemical modifications (m6A, pseudouridine, inosine) in the raw signal, because the native RNA molecule itself passes through the pore. m6A hypomethylation of specific transcripts has been linked to acute myeloid leukemia aggressiveness; detecting it directly from patient RNA samples could guide stratification. Short-read RNA-seq converts RNA to cDNA, erasing all modifications.

In all three cases, the extra information from long, direct-RNA reads provides therapeutically actionable detail that short reads cannot capture.
