INFO 636 - Assignment 1: De novo assembly of Escherichia coli with Velvet and Oases
Author: Vignu Gutta (vgutta@iu.edu)
System: IU Big Red 200 (Slurm, account c02374); work directory /N/scratch/vgutta/assembly
Data:   SRA run SRR21904868 (Illumina paired-end, 2 x 151 bp, 2,082,060 read pairs)

==========================================================================================
OVERVIEW
==========================================================================================
Goal: reconstruct the E. coli genome (~4.6-4.8 Mb) from short Illumina reads without a
reference genome (de novo assembly) with two de Bruijn graph tools:
  - Velvet 1.2.10 : genome assembler            -> contigs.fa
  - Oases 0.2.09  : transcriptome assembler built on Velvet's graph -> transcripts.fa
Both tools were run at 13 k-mer sizes (k = 21..141, step 10) as a Slurm job array, so all
assemblies ran in parallel. Velvet was also rerun at k = 21..71 with the expected coverage
set by hand. Assemblies were compared with a custom statistics script and with QUAST.

Result: Velvet k = 91 is optimal (215 contigs, 4.78 Mb, N50 118,786 bp, longest 397,292 bp).

Files in this submission:
  README.txt              this file
  assemble_array.sh       Slurm array: one task per (tool, k), 26 tasks
  summarize.sh            Slurm job run after the array: summary table, best k, contig lengths
  submit_all.sh           submits assemble_array.sh, then summarize.sh with a dependency
  assemble_manualcov.sh   Slurm array: Velvet at k = 21..71 with manual exp_cov / cov_cutoff
  assembly_stats.py       assembly statistics (N50, L50, ...) and per-contig lengths
  results/                summary.tsv, read_info.tsv, best_k.tsv, velvet_k91_lengths.tsv,
                          oases_k81_lengths.tsv, velvet Log files
  quast_results/          QUAST reports for all Velvet and all Oases assemblies
  report                  written report with results and discussion

==========================================================================================
STEP 0. ENVIRONMENT SETUP (Big Red 200 login node)
==========================================================================================
  module load conda                                   # conda/26.3.2
  conda create -n info636 -c conda-forge -c bioconda sra-tools velvet oases -y
  conda activate info636
  conda install -n info636 -c bioconda -c conda-forge quast -y
  pip install setuptools        # QUAST 5.3.0 imports distutils, removed in Python >= 3.12;
                                # setuptools provides a compatible distutils

Versions: Velvet 1.2.10 (MAXKMERLENGTH=191, OPENMP), Oases 0.2.09, sra-tools 3.4.1,
          QUAST 5.3.0

  mkdir -p /N/scratch/vgutta/assembly && cd /N/scratch/vgutta/assembly

==========================================================================================
STEP 1. DOWNLOAD THE READS (SRA Toolkit)
==========================================================================================
  mkdir -p ecoli_assembly/reads && cd ecoli_assembly/reads
  prefetch SRR21904868
  fasterq-dump --split-3 SRR21904868
  -> SRR21904868_1.fastq, SRR21904868_2.fastq (759 MB each), 4,164,120 reads

==========================================================================================
STEP 2. INSPECT THE READS
==========================================================================================
  head -4 SRR21904868_1.fastq       # read length 151; Nextera adapter read-through seen
  (read statistics are also computed by summarize.sh -> results/read_info.tsv)
  total bases 628,782,120 -> ~137x coverage of a 4.6 Mb genome
  No trimming was done; reads were assembled as downloaded.

==========================================================================================
STEP 3. ASSEMBLE: VELVET AND OASES AT 13 k VALUES (job array)
==========================================================================================
  cd /N/scratch/vgutta/assembly
  bash submit_all.sh
    -> sbatch --parsable assemble_array.sh                (26 tasks, 4 cores, 16 GB, 1 h each)
    -> sbatch --dependency=afterany:<array id> summarize.sh

k values: 21 31 41 51 61 71 81 91 101 111 121 131 141
Tasks 0-12 = Velvet, tasks 13-25 = Oases, same k list.

Velvet, for each k:
  velveth velvet/k<K> <K> -shortPaired -fastq -separate reads/SRR21904868_1.fastq reads/SRR21904868_2.fastq
  velvetg velvet/k<K> -exp_cov auto -cov_cutoff auto

Oases, for each k:
  velveth oases/k<K> <K> -shortPaired -fastq -separate reads/SRR21904868_1.fastq reads/SRR21904868_2.fastq
  velvetg oases/k<K> -read_trkg yes
  oases   oases/k<K>

Parameters:
  <K>                 k-mer (hash) length, odd, < read length (151), <= 191
  -shortPaired        short paired-end library
  -fastq -separate    FASTQ input, mates in two separate files
  -exp_cov auto       velvet estimates expected k-mer coverage (used for repeat
                      resolution and paired-end scaffolding)
  -cov_cutoff auto    removes nodes with coverage below half the estimated coverage
  -read_trkg yes      velvetg tracks reads per node; required by Oases
  oases               default settings
Each command was wrapped in /usr/bin/time -f "%e %M" to record wall time and peak memory.
k = 21 did not produce an assembly for either tool (see report).

==========================================================================================
STEP 4. VELVET WITH MANUAL COVERAGE (k = 21..71)
==========================================================================================
Velvet's automatic estimate failed at small k (k = 51: estimated coverage 1.04, cutoff
0.52, 1,864,632 nodes, only 5% of reads used). Rerun with coverage scaled from the
k = 91 estimate (26.6):
  exp_cov(k)  = 26.6 * (151 - k + 1) / (151 - 91 + 1)
  cov_cutoff  = exp_cov / 2
  sbatch assemble_manualcov.sh
    velveth velvet_manualcov/k<K> <K> -shortPaired -fastq -separate R1 R2
    velvetg velvet_manualcov/k<K> -exp_cov <exp_cov> -cov_cutoff <cov_cutoff>
  e.g. k = 61: -exp_cov 39.7 -cov_cutoff 19.9;  k = 71: -exp_cov 35.3 -cov_cutoff 17.6
Results: results/rows/velvetmcov_k<K>.tsv

==========================================================================================
STEP 5. ASSEMBLY STATISTICS AND OPTIMAL k (summarize.sh, assembly_stats.py)
==========================================================================================
  python3 assembly_stats.py summary <tool> <K> <fasta> <wall_s> <max_rss_kb>
    columns: tool, k, num_seqs, total_bp, longest_bp, N50, L50, num_ge_1kb,
             num_loci (Oases), wall_s, peak_mem_MB
  N50 = length L such that sequences >= L hold half of the assembled bases
  Lengths are computed from the sequences, not the FASTA headers
  (Velvet's NODE_n_length_L header gives L in k-mers; true length = L + k - 1).

Optimal k = highest N50 with total length close to the genome size:
  results/best_k.tsv : velvet 91, oases 81
  python3 assembly_stats.py lengths velvet/k91/contigs.fa     > results/velvet_k91_lengths.tsv
  python3 assembly_stats.py lengths oases/k81/transcripts.fa  > results/oases_k81_lengths.tsv
Full table: results/summary.tsv (also printed in logs/summary_<jobid>.log)

==========================================================================================
STEP 6. QUAST EVALUATION (no reference genome)
==========================================================================================
  KS="31 41 51 61 71 81 91 101 111 121 131 141"
  V=""; O=""; L=""
  for k in $KS; do V="$V ecoli_assembly/velvet/k$k/contigs.fa"; O="$O ecoli_assembly/oases/k$k/transcripts.fa"; L="$L,k$k"; done
  L=${L#,}
  quast $V -l "$L" -o quast_results/velvet --min-contig 200 --est-ref-size 4600000 -t 2
  quast $O -l "$L" -o quast_results/oases  --min-contig 200 --est-ref-size 4600000 -t 2

  --min-contig 200        ignore sequences shorter than 200 bp
  --est-ref-size 4600000  expected genome size, so QUAST also reports NG50 / LG50

==========================================================================================
HOW TO REPRODUCE
==========================================================================================
1. Steps 0-1 on the login node.
2. Put all scripts in /N/scratch/vgutta/assembly; run: bash submit_all.sh
3. Monitor: squeue -u $USER ; ls ecoli_assembly/results/rows/ | wc -l
4. sbatch assemble_manualcov.sh
5. Step 6 (QUAST) once the assemblies exist.
6. Copy results to home before scratch is purged:
   cp -r ecoli_assembly/results quast_results ~/assignment1_results/
