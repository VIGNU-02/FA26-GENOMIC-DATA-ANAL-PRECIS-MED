#!/bin/bash
#SBATCH -J align
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-11
#SBATCH -o /N/scratch/vgutta/asg2/aligned/align_%A_%a.log
module load conda
conda activate asg2
IDS=(SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 SRR22269877 \
     SRR22269878 SRR22269879 SRR22269880 SRR22269881 SRR22269882 SRR22269883)
id=${IDS[$SLURM_ARRAY_TASK_ID]}
WORK=/N/scratch/vgutta/asg2
hisat2 -p 8 -x $WORK/reference/grch38_index \
  -U $WORK/trimmed/${id}_trimmed.fq \
  --no-spliced-alignment -k 10 \
  --summary-file $WORK/aligned/${id}_summary.txt \
  | samtools sort -@ 4 -o $WORK/aligned/${id}.bam -
samtools index $WORK/aligned/${id}.bam
echo "ALIGN DONE: $id"
