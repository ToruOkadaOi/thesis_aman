#!/bin/bash
#SBATCH --job-name=stripediff_chr
#SBATCH --output=logs/stripediff_chr_%A_%a.out
#SBATCH --error=logs/stripediff_chr_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=medium
#SBATCH --array=0-23

# Set up environment
source ~/.bashrc
conda activate stripediff_aman

# Define chromosomes
chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
chr=${chroms[$SLURM_ARRAY_TASK_ID]}

# Run stripeDiff
echo "Running stripeDiff for $chr"
sh /usr/users/papantonis1/stripeDiff-master/src/stripeDiff.sh \
  -a /usr/users/papantonis1/stripeDiff-master/ctrl_${chr}_5kb.txt \
  -b /usr/users/papantonis1/stripeDiff-master/rbp1_${chr}_5kb.txt \
  -n control,rbp1,${chr} \
  -o /usr/users/papantonis1/stripeDiff-master/output_${chr}
