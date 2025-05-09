#!/bin/bash
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --job-name=microc_map_dedup
#SBATCH -o logs/microc_map_dedup_%A_%a.out
#SBATCH -e logs/microc_map_dedup_%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --array=0-1

set -euo pipefail

source ~/.bashrc
conda activate /usr/users/papantonis1/anaconda3/envs/microc_env

samples=(CTRL RBP1)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

WD="/usr/users/papantonis1/aman/microc_project"
FASTQ="$WD/00_raw_fastq/trimmed"
REF="$WD/refgen/hg38.fa"
CHROMS="$WD/refgen/hg38.chrom.sizes"
OUTDIR="$WD/03_pairs_bam/$sample"
mkdir -p "$OUTDIR"

echo "[1/5] Mapping and parsing $sample"
bwa mem -5SP -T0 -t 16 "$REF" \
  "$FASTQ/${sample}_R1.trimmed.fastq.gz" \
  "$FASTQ/${sample}_R2.trimmed.fastq.gz" | \
pairtools parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --nproc-in 8 \
  --nproc-out 8 \
  --chroms-path "$CHROMS" \
  - > "$OUTDIR/${sample}.parsed.pairsam"

echo "[2/5] Sorting $sample"
pairtools sort \
  --nproc 16 \
  --tmpdir "$OUTDIR/tmp" \
  -o "$OUTDIR/${sample}.sorted.pairsam" \
  "$OUTDIR/${sample}.parsed.pairsam"

rm -f "$OUTDIR/${sample}.parsed.pairsam"

echo "[3/5] Deduplicating $sample"
pairtools dedup \
  --nproc 16 \
  --output-stats "$OUTDIR/${sample}.dedup.stats" \
  -o "$OUTDIR/${sample}.dedup.pairsam" \
  "$OUTDIR/${sample}.sorted.pairsam"

rm -f "$OUTDIR/${sample}.sorted.pairsam"

echo "[4/5] Splitting to .pairs.gz and .bam"
pairtools split \
  --nproc 16 \
  --output-pairs - \
  --output-sam - \
  "$OUTDIR/${sample}.dedup.pairsam" | \
tee >(bgzip > "$OUTDIR/${sample}.pairs.gz") | \
samtools view -@ 16 -Sb -o "$OUTDIR/${sample}.bam" -

rm -f "$OUTDIR/${sample}.dedup.pairsam"

echo "[5/5] Indexing BAM"
samtools index "$OUTDIR/${sample}.bam"

echo "✅ Done for $sample — BAM and .pairs.gz generated and indexed."

