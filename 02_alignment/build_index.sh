#!/bin/bash
#SBATCH -J hisat_index
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH -o /N/scratch/vgutta/asg2/reference/index_%j.log
module load conda
conda activate asg2
cd /N/scratch/vgutta/asg2/reference
hisat2-build -p 16 Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38_index
echo "INDEX DONE"
