#!/bin/bash
#SBATCH -J count_ens
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o /N/scratch/vgutta/asg2/counts/count_ens_%j.log
module load conda
conda activate asg2
WORK=/N/scratch/vgutta/asg2
featureCounts -T 8 -M -O \
  -a $WORK/reference/Homo_sapiens.GRCh38.115.gtf \
  -t exon -g gene_id \
  -o $WORK/counts/ensembl_gene_counts.txt \
  $WORK/aligned/*.bam
echo "ENSEMBL COUNT DONE"
