INFO-B636 Assignment 1 - De novo assembly of E. coli with Velvet and Oases
Vignesh Gutta (vgutta@iu.edu)
Ran on IU Big Red 200, Slurm account c02374, working folder /N/scratch/vgutta/assembly
Data: SRR21904868, Illumina paired-end 2 x 151 bp, 2,082,060 read pairs

Summary
Velvet and Oases were each run at k = 21, 31, ..., 141. Velvet was also rerun at
k = 21-71 with coverage set by hand. Best assembly: Velvet k = 91 (215 contigs,
4.78 Mb, N50 118,786 bp, longest contig 397,292 bp).

Files
README.md               report with all results
README.txt              this file (steps and commands)
assemble_array.sh       Slurm array job, one task per tool and k (26 tasks)
summarize.sh            runs after the array, makes the summary table and contig lists
submit_all.sh           submits the array and then summarize.sh
assemble_manualcov.sh   Velvet at k = 21-71 with manual exp_cov / cov_cutoff
assembly_stats.py       N50, L50, totals and contig lengths
summary.tsv             stats for every assembly
velvet_k91_lengths.tsv  contig lengths, best Velvet assembly
oases_k81_lengths.tsv   transcript lengths, best Oases assembly
velvet_summary/      QUAST output for the Velvet assemblies
oases_summary/       QUAST output for the Oases assemblies


1) Setup (login node)

module load conda
conda create -n info636 -c conda-forge -c bioconda sra-tools velvet oases -y
conda activate info636
conda install -n info636 -c bioconda -c conda-forge quast -y
pip install setuptools     # quast 5.3.0 needs distutils, which python 3.12+ dropped

Versions: Velvet 1.2.10 (MAXKMERLENGTH=191, OpenMP), Oases 0.2.09,
sra-tools 3.4.1, QUAST 5.3.0


2) Download reads

mkdir -p /N/scratch/vgutta/assembly/ecoli_assembly/reads
cd /N/scratch/vgutta/assembly/ecoli_assembly/reads
prefetch SRR21904868
fasterq-dump --split-3 SRR21904868
head -4 SRR21904868_1.fastq

This gives SRR21904868_1.fastq and _2.fastq (4,164,120 reads, 151 bp).
Total is about 629 Mb, so roughly 137x coverage of a 4.6 Mb genome.
Some reads contain Nextera adapter sequence. I did not trim the reads.


3) Velvet and Oases for each k (job array)

cd /N/scratch/vgutta/assembly
bash submit_all.sh

k values: 21 31 41 51 61 71 81 91 101 111 121 131 141
Tasks 0-12 run Velvet, tasks 13-25 run Oases. Each task: 4 cores, 16 GB, 1 hour.

Velvet:
  velveth velvet/k<K> <K> -shortPaired -fastq -separate R1.fastq R2.fastq
  velvetg velvet/k<K> -exp_cov auto -cov_cutoff auto

Oases:
  velveth oases/k<K> <K> -shortPaired -fastq -separate R1.fastq R2.fastq
  velvetg oases/k<K> -read_trkg yes
  oases   oases/k<K>

Parameters
  <K>                 k-mer length (odd, less than the 151 bp read length)
  -shortPaired        paired-end short reads
  -fastq -separate    fastq input with the two mates in separate files
  -exp_cov auto       velvet estimates the expected k-mer coverage
  -cov_cutoff auto    drops nodes below half of that coverage (errors)
  -read_trkg yes      keeps track of reads per node, Oases needs this

Every command was run through /usr/bin/time to get wall time and peak memory.
k = 21 ran out of memory (16 GB) for both tools after about 11 minutes.


4) Velvet with coverage set by hand

At small k velvet's auto estimate was wrong (k = 51 estimated coverage 1.04,
so almost nothing was filtered). I reran k = 21-71 with:

  exp_cov    = 26.6 * (151 - k + 1) / (151 - 91 + 1)
  cov_cutoff = exp_cov / 2

26.6 is velvet's own estimate at k = 91.

sbatch assemble_manualcov.sh
  velvetg velvet_manualcov/k<K> -exp_cov <value> -cov_cutoff <value>
  e.g. k = 61: -exp_cov 39.7 -cov_cutoff 19.9


5) Stats and choosing the best k

python3 assembly_stats.py summary <tool> <K> <fasta> <wall_s> <max_rss_kb>
python3 assembly_stats.py lengths velvet/k91/contigs.fa > velvet_k91_lengths.tsv
python3 assembly_stats.py lengths oases/k81/transcripts.fa > oases_k81_lengths.tsv

Contig lengths are counted from the sequences. Velvet's NODE_x_length_L header
gives L in k-mers, so the real length is L + k - 1.
Best k = highest N50 with total length still close to the genome size:
Velvet k = 91, Oases k = 81.


6) QUAST

KS="31 41 51 61 71 81 91 101 111 121 131 141"
V=""; O=""; L=""
for k in $KS; do
  V="$V ecoli_assembly/velvet/k$k/contigs.fa"
  O="$O ecoli_assembly/oases/k$k/transcripts.fa"
  L="$L,k$k"
done
L=${L#,}
quast $V -l "$L" -o quast_results/velvet --min-contig 200 --est-ref-size 4600000 -t 2
quast $O -l "$L" -o quast_results/oases --min-contig 200 --est-ref-size 4600000 -t 2

k = 21 is left out because those runs did not finish.
--est-ref-size 4600000 lets QUAST report NG50 without a reference genome.
