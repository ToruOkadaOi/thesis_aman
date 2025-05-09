#!/bin/bash
#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --job-name=fastp_trim
#SBATCH -o logs/trim_%A_%a.out
#SBATCH -e logs/trim_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --array=0-1

set -euo pipefail

# Activate conda
source ~/.bashrc
conda activate /usr/users/papantonis1/anaconda3/envs/microc_env

# Sample names
samples=(CTRL RBP1)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Paths
WD="/usr/users/papantonis1/aman/microc_project"
RAW="$WD/00_raw_fastq"
TRIMMED="$RAW/trimmed"
REPORTS="$WD/01_qc/fastp_trimmed_reports"
mkdir -p "$TRIMMED" "$REPORTS"

echo "Trimming adapters for $sample..."

fastp \
  -i "$RAW/${sample}_R1.fastq.gz" \
  -I "$RAW/${sample}_R2.fastq.gz" \
  -o "$TRIMMED/${sample}_R1.trimmed.fastq.gz" \
  -O "$TRIMMED/${sample}_R2.trimmed.fastq.gz" \
  --detect_adapter_for_pe \
  -j "$REPORTS/${sample}.json" \
  -h "$REPORTS/${sample}.html" \
  -w 8

echo "✅ Trimming complete for $sample"

