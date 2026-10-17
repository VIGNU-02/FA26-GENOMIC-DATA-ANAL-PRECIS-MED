# Assignment 1 | Genome Assembly Using Velvet and Oases

Programmer Name: Vignesh Gutta  
Language: Bash (Slurm) and Python 3  
Date: 10/16/2026

## Description

For this assignment I assembled the *E. coli* genome de novo from Illumina short reads (SRA run SRR21904868) with Velvet and Oases, and looked at how the k-mer size and Velvet's coverage settings change the result.

I ran both tools at k = 21 to 141 (step 10) on Big Red 200, using a Slurm job array so the 26 assemblies ran at the same time. I also reran Velvet at small k with the coverage set manually. To compare assemblies I used a small Python script (contigs, total length, longest contig, N50, L50, time, memory) and QUAST.

Velvet with k = 91 gave the best assembly: 215 contigs, 4.78 Mb, N50 118,786 bp, longest contig 397,292 bp. Oases gave about 2.7x more sequence than the genome actually has, so Velvet is the better tool for this data.

## Requirements

| Tool | Version | Purpose |
| --- | --- | --- |
| SRA Toolkit | 3.4.1 | download reads (`prefetch`, `fasterq-dump`) |
| Velvet | 1.2.10 (MAXKMERLENGTH = 191, OpenMP) | genome assembly |
| Oases | 0.2.09 | transcriptome assembly on Velvet's graph |
| QUAST | 5.3.0 | assembly evaluation |
| Python | 3 | `assembly_stats.py` (N50 and contig lengths) |
| Conda | 26.3.2 | environment management |
| Slurm | Big Red 200, account c02374 | job arrays |

## Files in this repository

| File | Purpose |
| --- | --- |
| `README.md` | this report: steps, commands, results |
| `assemble_array.sh` | Slurm array, one task per (tool, k), 26 tasks |
| `summarize.sh` | runs after the array: summary table, best k, contig lengths |
| `submit_all.sh` | submits the array, then `summarize.sh` with a dependency |
| `assemble_manualcov.sh` | Slurm array: Velvet at k = 21-71 with manual coverage |
| `assembly_stats.py` | computes N50, L50, totals and every contig length |
| `summary.tsv` | statistics for every assembly |
| `velvet_k91_lengths.tsv` | every contig length, optimal Velvet assembly |
| `oases_k81_lengths.tsv` | every transcript length, optimal Oases assembly |
| `velvet_summary.zip` | QUAST report for all Velvet assemblies |
| `oases_summary.zip` | QUAST report for all Oases assemblies |

## Steps to Execute the Project

### Step 1: Environment setup (Big Red 200 login node)

```bash
module load conda
conda create -n info636 -c conda-forge -c bioconda sra-tools velvet oases -y
conda activate info636
conda install -n info636 -c bioconda -c conda-forge quast -y
pip install setuptools   # QUAST 5.3.0 imports distutils, removed in Python >= 3.12
```

### Step 2: Working directory

```bash
mkdir -p /N/scratch/vgutta/assembly/ecoli_assembly/reads
cd /N/scratch/vgutta/assembly
```

### Step 3: Download the reads

```bash
cd ecoli_assembly/reads
prefetch SRR21904868
fasterq-dump --split-3 SRR21904868      # -> SRR21904868_1.fastq, SRR21904868_2.fastq
head -4 SRR21904868_1.fastq             # 151 bp reads
cd ../..
```

### Step 4: Velvet and Oases at 13 k values (parallel job array)

```bash
bash submit_all.sh
# = sbatch --parsable assemble_array.sh          (26 tasks: 4 cores, 16 GB, 1 h each)
#   sbatch --dependency=afterany:<array id> summarize.sh
```

Commands run by each array task (`R1`, `R2` = the two FASTQ files):

```bash
# Velvet (tasks 0-12)
velveth velvet/k<K> <K> -shortPaired -fastq -separate R1 R2
velvetg velvet/k<K> -exp_cov auto -cov_cutoff auto

# Oases (tasks 13-25)
velveth oases/k<K> <K> -shortPaired -fastq -separate R1 R2
velvetg oases/k<K> -read_trkg yes
oases   oases/k<K>
```

| Parameter | Meaning |
| --- | --- |
| `<K>` | k-mer length: 21, 31, ..., 141 (odd, below the 151 bp read length) |
| `-shortPaired -fastq -separate` | paired-end FASTQ, mates in two files |
| `-exp_cov auto` | Velvet estimates expected k-mer coverage (repeat resolution, scaffolding) |
| `-cov_cutoff auto` | removes nodes below half the estimated coverage (errors) |
| `-read_trkg yes` | tracks reads per graph node; required by Oases |

Each command was wrapped in `/usr/bin/time` to record wall time and peak memory.

### Step 5: Velvet with manual coverage (k = 21-71)

```bash
sbatch assemble_manualcov.sh
# velvetg velvet_manualcov/k<K> -exp_cov <E> -cov_cutoff <E/2>
# E = 26.6 * (151 - k + 1) / (151 - 91 + 1)   (scaled from Velvet's own k = 91 estimate)
```

### Step 6: Statistics and optimal k

```bash
python3 assembly_stats.py summary <tool> <K> <fasta> <wall_s> <max_rss_kb>
python3 assembly_stats.py lengths velvet/k91/contigs.fa    > velvet_k91_lengths.tsv
python3 assembly_stats.py lengths oases/k81/transcripts.fa > oases_k81_lengths.tsv
```

Lengths are computed from the sequences: Velvet's `NODE_n_length_L` header gives L in k-mers, and the true length is L + k - 1. The optimal k is the highest N50 whose total length stays close to the genome size.

### Step 7: QUAST

```bash
KS="31 41 51 61 71 81 91 101 111 121 131 141"
V=""; O=""; L=""
for k in $KS; do V="$V ecoli_assembly/velvet/k$k/contigs.fa"; O="$O ecoli_assembly/oases/k$k/transcripts.fa"; L="$L,k$k"; done
L=${L#,}
quast $V -l "$L" -o quast_results/velvet --min-contig 200 --est-ref-size 4600000 -t 2
quast $O -l "$L" -o quast_results/oases  --min-contig 200 --est-ref-size 4600000 -t 2
```

## Results

### Data

| Property | Value |
| --- | --- |
| Layout | paired-end, 2 x 151 bp |
| Reads | 4,164,120 (2,082,060 pairs) |
| Total bases | 628,782,120 |
| Coverage of a 4.6 Mb genome | about 137x |

Some reads run into the Nextera adapter (`CTGTCTCTTATACACATCT`) and base quality drops toward the 3' end. Reads were assembled untrimmed.

### Velvet k sweep (automatic coverage)

| k | Contigs | Total (bp) | Longest (bp) | N50 (bp) | Time (s) | Peak memory (MB) |
| --- | --- | --- | --- | --- | --- | --- |
| 21 | out of memory (16 GB) | | | | 669 | 16,029 |
| 31 | 218,230 | 17,706,741 | 268 | 81 | 1,387 | 14,635 |
| 41 | 138,569 | 13,826,260 | 255 | 98 | 905 | 8,581 |
| 51 | 37,710 | 4,382,424 | 289 | 113 | 543 | 8,356 |
| 61 | 18,695 | 2,570,182 | 310 | 132 | 407 | 8,804 |
| 71 | 20,808 | 3,630,902 | 530 | 169 | 373 | 8,949 |
| 81 | 244 | 4,782,670 | 365,843 | 112,568 | 343 | 8,774 |
| 91 | 215 | 4,783,917 | 397,292 | 118,786 | 223 | 8,327 |
| 101 | 283 | 4,801,948 | 311,964 | 101,557 | 194 | 7,613 |
| 111 | 329 | 4,806,749 | 166,879 | 65,065 | 188 | 6,679 |
| 121 | 743 | 4,831,581 | 62,325 | 14,117 | 136 | 5,530 |
| 131 | 3,492 | 4,159,776 | 10,090 | 1,561 | 113 | 4,177 |
| 141 | 17 | 38,117 | 10,368 | 3,659 | 85 | 2,612 |

Quality rises to a peak at k = 91 and falls after it. At k = 131-141 k is too close to the read length: a 151 bp read holds only 151 - k + 1 k-mers (11 at k = 141), so one error removes most of them. At k = 21 both tools exceeded the 16 GB memory limit after about 11 minutes, because smaller k gives a much larger graph.

### Why small k failed: Velvet's coverage estimate

| Run | Estimated coverage | Cutoff | Graph nodes | Reads used |
| --- | --- | --- | --- | --- |
| k = 51, auto | 1.04 | 0.52 | 1,864,632 | 208,493 (5%) |
| k = 91, auto | 26.56 | 13.28 | 353 | 2,899,695 (70%) |

At k = 51 Velvet estimated coverage at 1x instead of roughly 44x. Sequencing errors and adapter k-mers create so many tiny low-coverage nodes that the estimate followed them, so nothing was filtered. Setting coverage by hand fixed it:

| k | exp_cov | cov_cutoff | Contigs | Total (bp) | Longest (bp) | N50 manual (bp) | N50 auto (bp) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 41 | 48.4 | 24.2 | 522 | 4,754,837 | 177,464 | 36,270 | 98 |
| 51 | 44.0 | 22.0 | 390 | 4,754,982 | 229,512 | 51,750 | 113 |
| 61 | 39.7 | 19.9 | 414 | 4,752,775 | 178,365 | 44,706 | 132 |
| 71 | 35.3 | 17.6 | 455 | 4,753,900 | 178,280 | 34,928 | 169 |

N50 rose by up to 450 times, so the small-k failure was a parameter problem, not a k problem. Even so, no run beat k = 91, because longer k-mers resolve more repeats.

### Oases k sweep

| k | Transcripts | Loci | Total (bp) | Longest (bp) | N50 (bp) | Time (s) | Peak memory (MB) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 21 | out of memory (16 GB) | | | | | 677 | 15,877 |
| 31 | 5,915 | 412 | 23,486,617 | 89,286 | 15,542 | 1,231 | 10,642 |
| 41 | 16,132 | 102 | 14,566,468 | 454,206 | 9,529 | 942 | 7,595 |
| 51 | 1,524 | 143 | 18,130,793 | 183,243 | 44,772 | 632 | 8,356 |
| 61 | 24,113 | 1,422 | 16,691,773 | 849,797 | 2,413 | 561 | 8,805 |
| 71 | 2,378 | 133 | 11,346,768 | 540,708 | 25,582 | 513 | 8,953 |
| 81 | 662 | 9 | 12,560,123 | 1,116,670 | 255,449 | 289 | 8,770 |
| 91 | 505 | 13 | 10,994,680 | 746,412 | 131,682 | 253 | 8,324 |
| 101 | 778 | 42 | 11,634,789 | 826,841 | 118,807 | 217 | 7,610 |
| 111 | 741 | 226 | 7,100,483 | 153,347 | 35,401 | 178 | 6,677 |
| 121 | 1,640 | 1,452 | 5,517,954 | 30,787 | 5,668 | 150 | 5,531 |
| 131 | 5,719 | 5,633 | 4,523,279 | 10,489 | 846 | 130 | 4,180 |
| 141 | 86 | 85 | 70,349 | 10,368 | 767 | 90 | 2,615 |

Oases' best N50 is at k = 81, but that assembly is 12.56 Mb, 2.7 times the genome: 662 transcripts in only 9 loci.

### QUAST (contigs >= 200 bp, estimated genome size 4.6 Mb)

| Assembly | Contigs | Total (bp) | Largest (bp) | N50 | NG50 | GC (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Velvet k = 81 | 211 | 4,776,832 | 365,843 | 112,568 | 112,568 | 50.83 |
| Velvet k = 91 | 197 | 4,780,511 | 397,292 | 118,786 | 119,969 | 50.83 |
| Velvet k = 101 | 283 | 4,801,948 | 311,964 | 101,557 | 105,796 | 50.84 |
| Oases k = 81 | 498 | 12,534,410 | 1,116,670 | 255,449 | 920,170 | 50.65 |
| Oases k = 91 | 395 | 10,978,673 | 746,412 | 131,682 | 611,191 | 50.98 |

GC content of 50.6-51.0% matches *E. coli* (about 50.8%). NG50 does not correct Oases' redundancy; it inflates it, because 12.5 Mb of duplicated sequence reaches half of 4.6 Mb within a few transcripts. Total length shows the redundancy directly.

### Optimal assemblies: number and length of contigs

Velvet, k = 91: 215 contigs, 4,783,917 bp. 16 contigs of at least 100 kb hold 56% of the assembly, 29 of at least 50 kb hold 75%, and 100 of at least 1 kb hold 99%. Every contig length (bp), longest first:

| Ranks | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1-10 | 397,292 | 296,317 | 232,925 | 216,319 | 162,089 | 151,049 | 147,418 | 131,682 | 125,930 | 125,694 |
| 11-20 | 124,234 | 123,583 | 119,969 | 118,786 | 112,588 | 105,794 | 94,468 | 92,323 | 88,458 | 69,220 |
| 21-30 | 68,814 | 68,018 | 67,804 | 67,643 | 59,053 | 57,697 | 52,208 | 51,224 | 50,535 | 49,385 |
| 31-40 | 49,282 | 46,187 | 45,156 | 43,343 | 43,085 | 41,349 | 40,300 | 40,216 | 38,458 | 38,219 |
| 41-50 | 37,172 | 36,208 | 35,150 | 34,863 | 31,090 | 31,026 | 29,094 | 28,963 | 28,022 | 26,914 |
| 51-60 | 26,361 | 25,993 | 23,454 | 19,653 | 19,384 | 16,572 | 16,447 | 16,045 | 15,963 | 15,027 |
| 61-70 | 12,058 | 11,420 | 10,943 | 10,383 | 10,226 | 9,783 | 9,486 | 9,332 | 8,168 | 7,892 |
| 71-80 | 6,662 | 5,769 | 4,804 | 4,496 | 3,961 | 3,566 | 3,397 | 3,036 | 3,036 | 2,581 |
| 81-90 | 2,533 | 2,272 | 2,146 | 2,066 | 2,013 | 1,994 | 1,676 | 1,546 | 1,444 | 1,384 |
| 91-100 | 1,343 | 1,318 | 1,287 | 1,161 | 1,111 | 1,044 | 1,040 | 1,031 | 1,030 | 1,020 |
| 101-110 | 991 | 956 | 948 | 893 | 881 | 868 | 868 | 862 | 855 | 826 |
| 111-120 | 823 | 823 | 823 | 820 | 814 | 747 | 735 | 728 | 689 | 677 |
| 121-130 | 608 | 607 | 587 | 586 | 518 | 516 | 515 | 514 | 514 | 494 |
| 131-140 | 455 | 454 | 452 | 430 | 423 | 407 | 399 | 398 | 398 | 396 |
| 141-150 | 394 | 391 | 391 | 390 | 383 | 379 | 344 | 342 | 339 | 337 |
| 151-160 | 328 | 318 | 318 | 316 | 305 | 304 | 291 | 290 | 281 | 258 |
| 161-170 | 258 | 252 | 252 | 248 | 247 | 247 | 246 | 246 | 242 | 239 |
| 171-180 | 238 | 238 | 238 | 238 | 238 | 237 | 237 | 232 | 228 | 225 |
| 181-190 | 225 | 225 | 225 | 224 | 223 | 223 | 221 | 221 | 221 | 218 |
| 191-200 | 217 | 217 | 217 | 214 | 214 | 209 | 201 | 197 | 195 | 195 |
| 201-210 | 192 | 192 | 192 | 192 | 191 | 191 | 191 | 189 | 188 | 187 |
| 211-215 | 185 | 184 | 183 | 181 | 181 | | | | | |

Oases, k = 81: 662 transcripts in 9 loci, 12,560,123 bp. 209 are at least 1 kb, 140 at least 10 kb and 46 at least 50 kb; the longest is 1,116,670 bp. All 662 lengths are in `oases_k81_lengths.tsv`.

## Summary of Findings

- Velvet with k = 91 and automatic coverage (215 contigs, 4.78 Mb, N50 118.8 kb). Contiguity rises from k = 81 to a peak at 91 and declines from 101 upward; k = 131-141 collapses because k nears the 151 bp read length.
- Velvet's automatic estimate failed at k <= 71 (1x instead of ~44x at k = 51). Setting it by hand raised N50 at k = 51 from 113 bp to 51.8 kb, but k = 91 stayed best.
- Velvet is clearly the better choice here. Velvet reconstructs the genome once at the right size. Oases, built for transcriptomes, treats graph branches as splice isoforms and outputs the same regions repeatedly: 12.56 Mb (2.7x the genome) in 662 transcripts from 9 loci, and it was slower (289 s vs 223 s at each tool's best k). Its larger N50 reflects duplication, not a better genome.
- run time and memory fall as k rises; k = 21 needed more than 16 GB for both tools.
- GC content 50.6-51.0% matches *E. coli*.
