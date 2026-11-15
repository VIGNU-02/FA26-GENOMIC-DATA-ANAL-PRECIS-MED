ASSIGNMENT 2 - miRNA EXPRESSION ANALYSIS IN SARS-CoV-2 INFECTED CELLS

Programmer name: Vignesh Gutta
Course: INFO-B636, Genomic Data Analysis for Precision Medicine
Date: 11/15/2026

Language of the scripts: Unix / Bash (download, QC, trimming, alignment, counting)
and R (DESeq2 differential expression, GO enrichment)

System: IU Big Red 200, Slurm scheduler, account c02374, partition general
Working directory: /N/scratch/vgutta/asg2
Human genome build: GRCh38, Ensembl release 115
miRNA annotation: miRBase v23.1


OVERVIEW

Twelve small-RNA-seq libraries from cultivated human respiratory cells
(BioProject PRJNA901149) were analyzed to find microRNAs that change on
SARS-CoV-2 infection, and to find what those microRNAs regulate. Three
biological replicates per group, in a 2 x 2 design:

                  24 h                              72 h
  Mock        SRR22269881, 882, 883             SRR22269875, 876, 877
  SARS-CoV-2  SRR22269878, 879, 880             SRR22269872, 873, 874

The pipeline downloads the reads, checks quality, trims adapter, aligns to
GRCh38, counts reads over mature miRNA features, then tests for differential
expression and runs GO enrichment on the target genes of the miRNAs that moved.
Two comparisons were run: mock against infected across both time points, and
24 h against 72 h within the infected samples.


PART 1 - SOFTWARE AND ENVIRONMENT

  module load conda
  conda create -n asg2 -c bioconda -c conda-forge \
    sra-tools fastqc trim-galore hisat2 samtools subread
  conda activate asg2

R packages were installed into the same environment with BiocManager::install():
DESeq2, clusterProfiler, org.Hs.eg.db and multiMiR, plus ggplot2 and ggrepel
for the figures.

Versions as installed:

  conda                                 26.3.2
  sra-tools (prefetch, fasterq-dump)    3.4.1
  Trim Galore                           v2.5.0
  featureCounts (subread)               v2.1.1
  clusterProfiler                       v4.18.4
  multiMiR                              miRTarBase through multiMiR database
                                        2.4.0, updated 2024-08-28
  FastQC                                conda build in env asg2
  HISAT2                                conda build in env asg2
  samtools                              conda build in env asg2
  R                                     run inside env asg2

HISAT2, samtools, FastQC and R version strings were not captured while the jobs
were running, so they are left unstated rather than guessed. Everything else is
the version that was actually installed.

A separate conda environment pins these versions, so a later change to the
system modules cannot alter the result.


PART 2 - DIRECTORY STRUCTURE

  mkdir -p /N/scratch/vgutta/asg2/{sra,fastqc,trimmed,reference,aligned,counts}
  cd /N/scratch/vgutta/asg2

Keeping each stage in its own folder means the array jobs in later parts can
write their output without colliding, and a failed step can be rerun without
touching anything upstream of it.


PART 3 - DOWNLOAD THE RAW SEQUENCING DATA

  cd /N/scratch/vgutta/asg2/sra
  for id in SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 \
            SRR22269877 SRR22269878 SRR22269879 SRR22269880 SRR22269881 \
            SRR22269882 SRR22269883
  do
      prefetch $id
      fasterq-dump $id --outdir /N/scratch/vgutta/asg2/sra
  done

prefetch pulls the SRA archive, fasterq-dump converts it to FASTQ. All 12 runs
are single-end, 50 bp.


PART 4 - QUALITY CONTROL WITH FastQC

  fastqc --outdir /N/scratch/vgutta/asg2/fastqc /N/scratch/vgutta/asg2/sra/*.fastq

(the FastQC build in this environment takes --outdir, not -o)

Every raw read is exactly 50 bp with adapter starting partway through, which is
what a small-RNA library should look like: the insert is only about 22 nt, so
the rest of each read is adapter sequence. The per-base quality and duplication
warnings FastQC raises on these files are expected for small RNA and are not a
reason to discard anything.


PART 5 - FINDING THE RIGHT ADAPTER

This step is worth doing before trimming, because getting it wrong is silent.

The first trim run used Trim Galore's --small_rna preset, which assumes the
adapter TGGAATTCTCGG. Only 22.5% of reads matched. Counting both candidate
adapters directly in the raw reads showed why:

  grep -c "AGATCGGAAGAGC" SRR22269872.fastq      # 3,879,925 hits
  grep -c "TGGAATTCTCGG"  SRR22269872.fastq      # far fewer

These libraries carry the Illumina TruSeq small-RNA adapter AGATCGGAAGAGC, so
it is passed explicitly in Part 6. After that change, adapter detection rose to
94.7% and the residual count in the trimmed output dropped to zero.


PART 6 - TRIMMING WITH Trim Galore

Script: 01_trimming/trim_array.sh (Slurm array, 12 tasks, one per sample)

  sbatch 01_trimming/trim_array.sh

The command each task runs:

  trim_galore -q 20 --phred33 --length 16 \
    -a AGATCGGAAGAGC \
    /N/scratch/vgutta/asg2/sra/${id}.fastq \
    -o /N/scratch/vgutta/asg2/trimmed/

  -q 20        trim low-quality 3' ends
  --phred33    quality encoding
  --length 16  discard reads shorter than 16 nt after trimming
  -a           the adapter identified in Part 5

--length 16 is set low on purpose. A standard 20 bp minimum would throw away
real mature miRNAs, which run about 18 to 25 nt. After trimming, read lengths
peak at 22 nt, which is where mature miRNAs should sit.


PART 7 - DOWNLOAD THE GRCh38 REFERENCE (Ensembl)

  cd /N/scratch/vgutta/asg2/reference
  wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

The Ensembl gene GTF is also downloaded, but only for the annotation comparison
in Part 11:

  wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
  gunzip Homo_sapiens.GRCh38.115.gtf.gz


PART 8 - DOWNLOAD THE miRBASE ANNOTATION AND FIX THE CHROMOSOME NAMES

  wget https://www.mirbase.org/download/hsa.gff3

miRBase gives the coordinates of every mature miRNA and precursor hairpin on
GRCh38, which is what makes it possible to count reads per miRNA instead of
per gene.

One thing has to be fixed first. miRBase writes chromosomes as "chr1" while the
Ensembl genome writes them as "1", so the two files disagree on every
chromosome name. featureCounts does not warn about this; it simply matches
nothing and returns a table of zeros. Stripping the prefix fixes it:

  sed 's/^chr//' hsa.gff3 > hsa.ensembl.gff3

All counting from here uses hsa.ensembl.gff3.


PART 9 - BUILD THE HISAT2 GENOME INDEX

Script: 02_alignment/build_index.sh (16 cores, 64 GB, about an hour)

  sbatch 02_alignment/build_index.sh

  hisat2-build -p 16 Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38_index

Indexing turns the FASTA into a structure HISAT2 can search quickly. It depends
only on the genome, not on the reads, so this job was submitted while trimming
was still running.


PART 10 - ALIGNMENT TO GRCh38

Script: 02_alignment/align_array.sh (Slurm array, 12 tasks)

  sbatch 02_alignment/align_array.sh

The command each task runs:

  hisat2 -p 8 -x $WORK/reference/grch38_index \
    -U $WORK/trimmed/${id}_trimmed.fq \
    --no-spliced-alignment -k 10 \
    --summary-file $WORK/aligned/${id}_summary.txt \
    | samtools sort -@ 4 -o $WORK/aligned/${id}.bam -
  samtools index $WORK/aligned/${id}.bam

Two settings here differ from an mRNA run:

  --no-spliced-alignment   A mature miRNA is one short piece of sequence with no
                           exon junctions in it. Leaving splicing on lets HISAT2
                           split a 22 nt read across two distant loci, which can
                           only produce false alignments.

  -k 10                    Report up to 10 alignments per read. Members of a
                           miRNA family share most of their sequence, so a read
                           can map equally well to several real loci, and
                           forcing a single best hit would discard that.

Piping HISAT2 straight into samtools sort writes a sorted BAM in one pass, with
no intermediate SAM file. A SAM for 12 samples would have cost a few hundred GB
of scratch for nothing.

Overall alignment rates (02_alignment/alignment_rates.txt):

  SRR22269872  54.35%        SRR22269878  62.12%
  SRR22269873  61.06%        SRR22269879  52.93%
  SRR22269874  67.51%        SRR22269880  54.24%
  SRR22269875  65.15%        SRR22269881  60.47%
  SRR22269876  67.65%        SRR22269882  44.67%
  SRR22269877  72.88%        SRR22269883  23.72%

Rates in the 50 to 70% range are normal for small-RNA libraries, which carry
plenty of non-miRNA material alongside the miRNAs: tRNA and rRNA fragments and
degradation products. SRR22269883 at 23.72% is an outlier and is flagged in the
report.


PART 11 - COUNTING READS PER miRNA

Script: 03_counting/count.sh

  sbatch 03_counting/count.sh

  featureCounts -T 8 -M -O \
    -a $WORK/reference/hsa.ensembl.gff3 \
    -F GTF -t miRNA -g Name \
    -o $WORK/counts/mirna_counts.txt \
    $WORK/aligned/*.bam

  -F GTF       read the miRBase gff3 with the GTF parser
  -t miRNA     count over mature miRNA features, not precursor hairpins
  -g Name      group by miRNA name (hsa-let-7a-5p), not by accession
  -M           count multi-mapping reads
  -O           count reads that overlap more than one feature

-M and -O are both needed because miRNA family members share sequence and some
mature miRNAs sit inside overlapping annotations. They do inflate counts for
those families, which the report states.

Clean the table for R. This drops the featureCounts header, keeps the ID column
plus the 12 count columns, and strips the full BAM paths out of the sample names
so they match the metadata file:

  cd /N/scratch/vgutta/asg2/counts
  cut -f1,7-18 mirna_counts.txt | grep -v "^#" \
    | sed '1s|/N/scratch/vgutta/asg2/aligned/||g; 1s|\.bam||g' \
    > mirna_counts_clean.tsv

Comparison against the full Ensembl gene annotation, run only to measure how
much of that output is actually miRNA:

  sbatch 03_counting/count_ensembl.sh

  featureCounts -T 8 -M -O \
    -a $WORK/reference/Homo_sapiens.GRCh38.115.gtf \
    -t exon -g gene_id \
    -o $WORK/counts/ensembl_gene_counts.txt \
    $WORK/aligned/*.bam

The Ensembl GTF counts 78,899 genes, of which 1,879 (2.4%) are miRNA biotype.
It assigns roughly three times more reads than miRBase, but 97.6% of the
features it measures are not miRNAs, so miRBase is what the analysis uses.


PART 12 - SAMPLE METADATA

03_counting/metadata.csv, built from assignment_2_info.csv:

  sample,condition,time
  SRR22269872,SARS,72h
  SRR22269873,SARS,72h
  SRR22269874,SARS,72h
  SRR22269875,Mock,72h
  SRR22269876,Mock,72h
  SRR22269877,Mock,72h
  SRR22269878,SARS,24h
  SRR22269879,SARS,24h
  SRR22269880,SARS,24h
  SRR22269881,Mock,24h
  SRR22269882,Mock,24h
  SRR22269883,Mock,24h

Row order matches the column order of mirna_counts_clean.tsv. DESeq2 does not
check this for you, so the script verifies it before building the dataset.


PART 13 - DIFFERENTIAL EXPRESSION WITH DESeq2

Script: 04_differential_expression/deseq2_analysis.R

  Rscript 04_differential_expression/deseq2_analysis.R

What the script does:

  - reads mirna_counts_clean.tsv and metadata.csv and checks the sample names
    line up in the same order
  - pre-filters:  keep <- rowSums(counts >= 10) >= 3
    This keeps a miRNA only if at least 3 samples have 10 or more reads for it,
    which leaves 527 of 2,783. The cut is there because a miRNA seen in one
    sample out of twelve carries no usable information and still costs a
    multiple-testing penalty for every other miRNA.
  - comparison 1, all 12 samples, design ~ time + condition
      contrast c("condition","SARS","Mock")
    Putting time in the design controls for the time point so the condition
    term measures infection on its own.
  - comparison 2, the 6 infected samples only, design ~ time
      contrast c("time","72h","24h")
  - significance: subset(res, padj < 0.05)
  - rlog transform for the PCA. vst() fails on this dataset because too few
    features have a mean normalized count above 5, which is the nsub error it
    raises.

Outputs in 04_differential_expression/:

  DEG_Mock_vs_SARS_full.csv             all 527 tested
  DEG_Mock_vs_SARS_significant.csv      28 at padj < 0.05 (11 up, 17 down)
  DEG_SARS_24h_vs_72h_full.csv          all 527 tested
  DEG_SARS_24h_vs_72h_significant.csv   49 at padj < 0.05 (22 up, 27 down)
  miRNA_normalized_expression.csv       normalized counts, all samples
  PCA_plot.pdf

Significance is padj < 0.05 with no fold-change filter. Adding |log2FC| > 1
drops comparison 1 from 28 to 4. That cutoff comes from mRNA work, and miRNAs
tune their targets rather than switch them, so it removes miRNAs that are
statistically solid: miR-877-5p, for instance, at log2FC 0.71 and padj 6.1e-05.
Fold change is reported in every table so the size of each change stays visible.

Script: 04_differential_expression/volcano_plots.R (ggplot2 and ggrepel)

  Rscript 04_differential_expression/volcano_plots.R

  volcano_Mock_vs_SARS.png
  volcano_SARS_24h_vs_72h.png


PART 14 - GO ENRICHMENT

Script: 05_go_enrichment/go_final.R

  Rscript 05_go_enrichment/go_final.R

miRNAs have almost no direct GO annotation of their own, so the enrichment runs
on their target genes instead. Up- and down-regulated miRNAs are handled
separately, because a miRNA going up suppresses its targets while one going
down releases them, and putting both in one list cancels the two effects out.

Target retrieval, then a filter down to strongly supported interactions:

  get_multimir(org = "hsa", mirna = <miRNA list>,
               table = "mirtarbase", summary = FALSE)

  grepl("Functional MTI", support_type) & !grepl("Weak", support_type)

An earlier attempt requested TarBase as well and it never returned. TarBase
also contributes lower-confidence high-throughput calls, so miRTarBase alone is
the more conservative source.

The background set matters as much as the query set. enrichGO defaults to the
whole genome, which is the wrong comparison here: miRNA targets are not a random
sample of genes, they skew toward heavily studied ones, so against a genome
background almost everything comes out enriched. The background used instead is
the union of strongly validated targets of all 527 expressed miRNAs:

  expressed <- rownames(cnt)[rowSums(cnt >= 10) >= 3]        # 527 miRNAs
  bg <- get_multimir(org = "hsa", mirna = expressed,
                     table = "mirtarbase", summary = FALSE)@data
  bg <- bg[grepl("Functional MTI", bg$support_type) &
           !grepl("Weak", bg$support_type), ]
  universe <- unique(bg$target_symbol)                       # 2,454 genes

  ego <- enrichGO(gene = tg, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

That change cut the up-regulated mock-vs-SARS set from 1,963 terms to 13.

Outputs in 05_go_enrichment/:

  GOfinal_Mock_vs_SARS_UP.csv             13 terms
  GOfinal_Mock_vs_SARS_DOWN.csv           749 terms
  GOfinal_SARS_24h_vs_72h_UP.csv          307 terms
  GOfinal_SARS_24h_vs_72h_DOWN.csv        573 terms


PART 15 - RESULTS AND INTERPRETATION

The numbers, the miRNA family patterns, the autophagy enrichment and the
quality limitations of this dataset are written up in REPORT.md.


FILES IN THIS REPOSITORY

Folders follow the pipeline in order. Each one holds the script for that step
and the output that step produced.

  README.txt                              this file
  README.md                               overview
  REPORT.md                               results and interpretation

  01_trimming/
    trim_array.sh                         adapter and quality trimming

  02_alignment/
    build_index.sh                        HISAT2 index build
    align_array.sh                        alignment to GRCh38
    alignment_rates.txt                   overall alignment rate per sample

  03_counting/
    count.sh                              featureCounts against miRBase
    count_ensembl.sh                      featureCounts against the Ensembl GTF
    metadata.csv                          sample sheet
    mirna_counts_clean.tsv                raw count matrix, 2,783 x 12
    mirna_counts.txt.summary              featureCounts assignment summary

  04_differential_expression/
    deseq2_analysis.R                     DESeq2, both comparisons, PCA
    volcano_plots.R                       volcano plots
    DEG_Mock_vs_SARS_full.csv
    DEG_Mock_vs_SARS_significant.csv
    DEG_SARS_24h_vs_72h_full.csv
    DEG_SARS_24h_vs_72h_significant.csv
    miRNA_normalized_expression.csv
    PCA_plot.pdf
    volcano_Mock_vs_SARS.png
    volcano_SARS_24h_vs_72h.png

  05_go_enrichment/
    go_final.R                            GO enrichment with corrected background
    GOfinal_Mock_vs_SARS_UP.csv
    GOfinal_Mock_vs_SARS_DOWN.csv
    GOfinal_SARS_24h_vs_72h_UP.csv
    GOfinal_SARS_24h_vs_72h_DOWN.csv

The SRA downloads, FastQC reports, trimmed FASTQ files, BAM files and the
GRCh38 index are not in the repository because of their size. The commands
above regenerate all of them.


RUN ORDER

  1.  download the SRA runs                                  Part 3
  2.  fastqc                                                 Part 4
  3.  sbatch 01_trimming/trim_array.sh                       Part 6
  4.  sbatch 02_alignment/build_index.sh                     Part 9
  5.  sbatch 02_alignment/align_array.sh                     Part 10
  6.  sbatch 03_counting/count.sh                            Part 11
  7.  clean the count table                                  Part 11
  8.  Rscript 04_differential_expression/deseq2_analysis.R   Part 13
  9.  Rscript 04_differential_expression/volcano_plots.R     Part 13
  10. Rscript 05_go_enrichment/go_final.R                    Part 14

Steps 3 and 4 are independent, so the index build ran alongside trimming.

The scripts use absolute paths under /N/scratch/vgutta/asg2 and can be
submitted from any working directory. Change that prefix to run them elsewhere.
