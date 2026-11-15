ASSIGNMENT 2 - REPRODUCIBILITY NOTES
miRNA differential expression, SARS-CoV-2 infected vs mock human respiratory cells
Vignesh Gutta | INFO-B636 Genomic Data Analysis for Precision Medicine | 11/15/2026

====================================================================
1. ENVIRONMENT
====================================================================

System        : IU Big Red 200 (Cray Shasta), Slurm scheduler
Slurm account : c02374
Partition     : general
Working dir   : /N/scratch/vgutta/asg2
Shell         : bash
Package mgr   : conda 26.3.2 (loaded with `module load conda`)
Conda env     : asg2
R             : run on Big Red inside the asg2 env

Software versions as installed:
  prefetch / fasterq-dump (sra-tools) 3.4.1
  FastQC                              (conda build in env asg2)
  Trim Galore                         v2.5.0
  cutadapt                            (dependency of Trim Galore)
  HISAT2                              (conda build in env asg2; version not recorded at run time)
  samtools                            (conda build in env asg2; version not recorded at run time)
  featureCounts (subread)             v2.1.1
  R / DESeq2
  clusterProfiler                     v4.18.4
  org.Hs.eg.db
  multiMiR                            (miRTarBase via multiMiR database 2.4.0, updated 2024-08-28)
  ggplot2, ggrepel

Environment setup:

  module load conda
  conda create -n asg2 -c bioconda -c conda-forge \
    sra-tools fastqc trim-galore hisat2 samtools subread
  conda activate asg2

R packages were installed inside the same environment with
BiocManager::install() for DESeq2, clusterProfiler, org.Hs.eg.db and multiMiR.

====================================================================
2. REFERENCE DATA
====================================================================

Human genome build : GRCh38
Source             : Ensembl release 115
File               : Homo_sapiens.GRCh38.dna.primary_assembly.fa

  cd /N/scratch/vgutta/asg2/reference
  wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Gene annotation (used only for the miRBase-vs-Ensembl comparison):

  wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
  gunzip Homo_sapiens.GRCh38.115.gtf.gz

miRNA annotation : miRBase v23.1, GRCh38 coordinates

  wget https://www.mirbase.org/download/hsa.gff3

IMPORTANT - chromosome naming. miRBase writes chromosomes as "chr1", the
Ensembl genome writes them as "1". Without fixing this, featureCounts
matches nothing and every count comes out zero:

  sed 's/^chr//' hsa.gff3 > hsa.ensembl.gff3

All downstream counting uses hsa.ensembl.gff3.

====================================================================
3. SEQUENCE DATA
====================================================================

BioProject : PRJNA901149
12 single-end small-RNA-seq runs, 50 bp:

  SRR22269872  SARS  72h
  SRR22269873  SARS  72h
  SRR22269874  SARS  72h
  SRR22269875  Mock  72h
  SRR22269876  Mock  72h
  SRR22269877  Mock  72h
  SRR22269878  SARS  24h
  SRR22269879  SARS  24h
  SRR22269880  SARS  24h
  SRR22269881  Mock  24h
  SRR22269882  Mock  24h
  SRR22269883  Mock  24h

Download:

  cd /N/scratch/vgutta/asg2/sra
  for id in SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 \
            SRR22269877 SRR22269878 SRR22269879 SRR22269880 SRR22269881 \
            SRR22269882 SRR22269883; do
      prefetch $id
      fasterq-dump $id --outdir /N/scratch/vgutta/asg2/sra
  done

====================================================================
4. QUALITY CONTROL
====================================================================

  fastqc --outdir /N/scratch/vgutta/asg2/fastqc /N/scratch/vgutta/asg2/sra/*.fastq

(the FastQC build in this env takes --outdir, not -o)

Every raw read was exactly 50 bp with adapter starting partway through,
which is the expected shape for a small-RNA library: the insert is only
about 22 nt, so the rest of the read is adapter.

Finding the adapter. The first trim run used Trim Galore's --small_rna
preset, which assumes TGGAATTCTCGG, and only 22.5% of reads matched.
Checking the reads directly:

  grep -c "AGATCGGAAGAGC" SRR22269872.fastq     # 3,879,925 hits
  grep -c "TGGAATTCTCGG"  SRR22269872.fastq     # far fewer

The real adapter is the Illumina TruSeq small-RNA adapter AGATCGGAAGAGC,
so it is passed explicitly in trim_array.sh. After that change, adapter
detection was 94.7% and the residual count in the trimmed output was zero.

====================================================================
5. TRIMMING
====================================================================

Script: trim_array.sh   (Slurm array, 12 tasks)

  sbatch trim_array.sh

Core command per sample:

  trim_galore -q 20 --phred33 --length 16 \
    -a AGATCGGAAGAGC \
    /N/scratch/vgutta/asg2/sra/${id}.fastq \
    -o /N/scratch/vgutta/asg2/trimmed/

  -q 20        trim low-quality 3' ends
  --length 16  discard reads shorter than 16 nt after trimming
  -a           adapter sequence (see section 4)

Post-trim read lengths peaked at 22 nt.

====================================================================
6. GENOME INDEX
====================================================================

Script: build_index.sh   (16 cores, 64 GB, ~1 h)

  sbatch build_index.sh

  hisat2-build -p 16 Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38_index

====================================================================
7. ALIGNMENT
====================================================================

Script: align_array.sh   (Slurm array, 12 tasks)

  sbatch align_array.sh

Core command per sample:

  hisat2 -p 8 -x $WORK/reference/grch38_index \
    -U $WORK/trimmed/${id}_trimmed.fq \
    --no-spliced-alignment -k 10 \
    --summary-file $WORK/aligned/${id}_summary.txt \
    | samtools sort -@ 4 -o $WORK/aligned/${id}.bam -
  samtools index $WORK/aligned/${id}.bam

  --no-spliced-alignment   mature miRNAs are single short pieces, there are
                           no exon junctions to find, and allowing splicing
                           only creates false alignments
  -k 10                    report up to 10 alignments per read; miRNA family
                           members share sequence, so a read can legitimately
                           map to more than one locus

Overall alignment rates (alignment_rates.txt):

  SRR22269872  54.35%      SRR22269878  62.12%
  SRR22269873  61.06%      SRR22269879  52.93%
  SRR22269874  67.51%      SRR22269880  54.24%
  SRR22269875  65.15%      SRR22269881  60.47%
  SRR22269876  67.65%      SRR22269882  44.67%
  SRR22269877  72.88%      SRR22269883  23.72%

====================================================================
8. COUNTING
====================================================================

Script: count.sh

  sbatch count.sh

  featureCounts -T 8 -M -O \
    -a $WORK/reference/hsa.ensembl.gff3 \
    -F GTF -t miRNA -g Name \
    -o $WORK/counts/mirna_counts.txt \
    $WORK/aligned/*.bam

  -F GTF       miRBase gff3 is read with the GTF parser
  -t miRNA     count over mature miRNA features, not precursor hairpins
  -g Name      group by miRNA name (hsa-let-7a-5p), not by accession
  -M           count multi-mapping reads
  -O           count reads overlapping more than one feature

-M and -O are both needed here because miRNA family members share
sequence and some mature miRNAs sit inside overlapping annotations.
They do inflate counts for those families, which is noted in the report.

Clean the count table for R (drop the featureCounts header, keep the ID
column plus the 12 count columns, strip full paths from the sample names):

  cd /N/scratch/vgutta/asg2/counts
  cut -f1,7-18 mirna_counts.txt | grep -v "^#" \
    | sed '1s|/N/scratch/vgutta/asg2/aligned/||g; 1s|\.bam||g' \
    > mirna_counts_clean.tsv

Comparison run against the full Ensembl gene annotation, used only to
quantify how much of that output is actually miRNA:

  sbatch count_ensembl.sh

  featureCounts -T 8 -M -O \
    -a $WORK/reference/Homo_sapiens.GRCh38.115.gtf \
    -t exon -g gene_id \
    -o $WORK/counts/ensembl_gene_counts.txt \
    $WORK/aligned/*.bam

Result: 78,899 genes counted, of which 1,879 (2.4%) are miRNA biotype.

====================================================================
9. SAMPLE METADATA
====================================================================

counts/metadata.csv, built from assignment_2_info.csv:

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

Row order matches the column order of mirna_counts_clean.tsv.

====================================================================
10. DIFFERENTIAL EXPRESSION
====================================================================

Script: deseq2_analysis.R

  Rscript deseq2_analysis.R

Steps:

  - read mirna_counts_clean.tsv and metadata.csv, check that sample
    names line up in the same order
  - pre-filter: keep <- rowSums(counts >= 10) >= 3        # 527 of 2,783 kept
  - comparison 1, all 12 samples, design ~ time + condition
      contrast c("condition","SARS","Mock")
  - comparison 2, 6 infected samples only, design ~ time
      contrast c("time","72h","24h")
  - significance: subset(res, padj < 0.05)
  - rlog transform for PCA (vst() fails on this dataset: too few features
    with mean normalized count > 5)

Outputs:
  DEG_Mock_vs_SARS_full.csv            all 527 tested
  DEG_Mock_vs_SARS_significant.csv     28 at padj < 0.05
  DEG_SARS_24h_vs_72h_full.csv         all 527 tested
  DEG_SARS_24h_vs_72h_significant.csv  49 at padj < 0.05
  miRNA_normalized_expression.csv      normalized counts, all samples
  PCA_plot.pdf

Script: volcano_plots.R   (ggplot2 + ggrepel)

  Rscript volcano_plots.R

Outputs:
  volcano_Mock_vs_SARS.png
  volcano_SARS_24h_vs_72h.png

====================================================================
11. GO ENRICHMENT
====================================================================

Script: go_final.R

  Rscript go_final.R

miRNAs have little direct GO annotation, so enrichment runs on their
validated target genes. Up- and down-regulated miRNAs are handled
separately, because a miRNA going up suppresses its targets while one
going down de-represses them.

Target retrieval:

  get_multimir(org = "hsa", mirna = <miRNA list>,
               table = "mirtarbase", summary = FALSE)

then kept only strongly supported interactions:

  grepl("Functional MTI", support_type) & !grepl("Weak", support_type)

(TarBase was requested in an earlier attempt and did not return; it also
contributes lower-confidence high-throughput calls, so miRTarBase alone
is the more conservative source.)

Background set. enrichGO defaults to the whole genome, which is wrong
here: miRNA targets are biased toward heavily studied genes, so almost
everything comes out enriched. The background used instead is the union
of strongly validated targets of all 527 expressed miRNAs:

  expressed <- rownames(cnt)[rowSums(cnt >= 10) >= 3]        # 527 miRNAs
  bg <- get_multimir(org = "hsa", mirna = expressed,
                     table = "mirtarbase", summary = FALSE)@data
  bg <- bg[grepl("Functional MTI", bg$support_type) &
           !grepl("Weak", bg$support_type), ]
  universe <- unique(bg$target_symbol)                       # 2,454 genes

Enrichment:

  ego <- enrichGO(gene = tg, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

Outputs:
  GOfinal_Mock_vs_SARS_UP.csv            13 terms
  GOfinal_Mock_vs_SARS_DOWN.csv          749 terms
  GOfinal_SARS_24h_vs_72h_UP.csv         307 terms
  GOfinal_SARS_24h_vs_72h_DOWN.csv       573 terms

====================================================================
12. FILES IN THIS REPOSITORY
====================================================================

Scripts
  trim_array.sh                        adapter and quality trimming (Slurm array)
  build_index.sh                       HISAT2 index build
  align_array.sh                       alignment to GRCh38 (Slurm array)
  count.sh                             featureCounts against miRBase
  count_ensembl.sh                     featureCounts against Ensembl GTF (comparison)
  deseq2_analysis.R                    DESeq2, both comparisons, PCA
  volcano_plots.R                      volcano plots
  go_final.R                           GO enrichment with corrected background

Data and results
  metadata.csv                         sample sheet
  mirna_counts_clean.tsv               raw count matrix, 2,783 x 12
  mirna_counts.txt.summary             featureCounts assignment summary
  miRNA_normalized_expression.csv      DESeq2 normalized counts
  alignment_rates.txt                  per-sample overall alignment rate
  DEG_Mock_vs_SARS_full.csv
  DEG_Mock_vs_SARS_significant.csv
  DEG_SARS_24h_vs_72h_full.csv
  DEG_SARS_24h_vs_72h_significant.csv
  GOfinal_Mock_vs_SARS_UP.csv
  GOfinal_Mock_vs_SARS_DOWN.csv
  GOfinal_SARS_24h_vs_72h_UP.csv
  GOfinal_SARS_24h_vs_72h_DOWN.csv

Figures
  PCA_plot.pdf
  volcano_Mock_vs_SARS.png
  volcano_SARS_24h_vs_72h.png

Write-up
  REPORT.md                            results and interpretation
  README.md                            overview
  README.txt                           this file

====================================================================
13. RUN ORDER
====================================================================

  1.  download SRA runs                   (section 3)
  2.  fastqc                              (section 4)
  3.  sbatch trim_array.sh                (section 5)
  4.  sbatch build_index.sh               (section 6, can run alongside step 3)
  5.  sbatch align_array.sh               (section 7)
  6.  sbatch count.sh                     (section 8)
  7.  clean the count table               (section 8)
  8.  Rscript deseq2_analysis.R           (section 10)
  9.  Rscript volcano_plots.R             (section 10)
  10. Rscript go_final.R                  (section 11)

Steps 3 and 4 are independent, so the index build was submitted while
trimming was still running.
