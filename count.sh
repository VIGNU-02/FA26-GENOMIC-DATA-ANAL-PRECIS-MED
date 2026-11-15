#!/bin/bash
#SBATCH -J count
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o /N/scratch/vgutta/asg2/counts/count_%j.log
module load conda
conda activate asg2
WORK=/N/scratch/vgutta/asg2
featureCounts -T 8 -M -O \
  -a $WORK/reference/hsa.ensembl.gff3 \
  -F GTF -t miRNA -g Name \
  -o $WORK/counts/mirna_counts.txt \
  $WORK/aligned/*.bam
echo "COUNT DONE"
