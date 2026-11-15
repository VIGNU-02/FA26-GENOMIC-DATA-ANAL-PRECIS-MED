#!/bin/bash
#SBATCH -J trim
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-11
#SBATCH -o /N/scratch/vgutta/asg2/trimmed/trim_%A_%a.log
module load conda
conda activate asg2
IDS=(SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 SRR22269877 \
     SRR22269878 SRR22269879 SRR22269880 SRR22269881 SRR22269882 SRR22269883)
id=${IDS[$SLURM_ARRAY_TASK_ID]}
trim_galore -q 20 --phred33 --length 16 \
  -a AGATCGGAAGAGC \
  /N/scratch/vgutta/asg2/sra/${id}.fastq \
  -o /N/scratch/vgutta/asg2/trimmed/
echo "TRIM DONE: $id"
