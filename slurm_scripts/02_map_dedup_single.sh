#!/bin/bash
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --job-name=microc
#SBATCH -o logs/microc_%x_%j.out
#SBATCH -e logs/microc_%x_%j.err
#SBATCH --time=46:00:00

set -euo pipefail

source ~/.bashrc
conda activate /usr/users/papantonis1/anaconda3/envs/microc_env

sample=$1

# Safety check
if [[ "$sample" != "CTRL" && "$sample" != "RBP1" ]]; then
  echo " ^}^l Invalid sample name: $sample"
  exit 1
fi

WD="/usr/users/papantonis1/aman/microc_project"
FASTQ="$WD/00_raw_fastq/trimmed"
REF="$WD/refgen/hg38.fa"
CHROMS="$WD/refgen/hg38.chrom.sizes"
OUTDIR="$WD/03_pairs_bam/$sample"
mkdir -p "$OUTDIR"

echo "[1/6] Mapping and parsing $sample"
bwa mem -5SP -T0 -t 8 "$REF" \
  "$FASTQ/${sample}_R1.trimmed.fastq.gz" \
  "$FASTQ/${sample}_R2.trimmed.fastq.gz" | \
pairtools parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --nproc-in 4 --nproc-out 4 \
  --chroms-path "$CHROMS" \
  - > "$OUTDIR/${sample}.parsed.pairsam"

echo "[2/6] Sorting $sample"
pairtools sort \
  --nproc 8 \
  --tmpdir "$OUTDIR/tmp" \
  -o "$OUTDIR/${sample}.sorted.pairsam" \
  "$OUTDIR/${sample}.parsed.pairsam"

rm "$OUTDIR/${sample}.parsed.pairsam"

echo "[3/6] Deduplicating $sample"
pairtools dedup \
  --nproc 8 \
  --output-stats "$OUTDIR/${sample}.dedup.stats" \
  -o "$OUTDIR/${sample}.dedup.pairsam" \
  "$OUTDIR/${sample}.sorted.pairsam"

rm "$OUTDIR/${sample}.sorted.pairsam"

echo "[4/6] Splitting to BAM and .pairs.gz"
pairtools split \
  --nproc 8 \
  --output-pairs - \
  --output-sam - \
  "$OUTDIR/${sample}.dedup.pairsam" | \
tee >(bgzip > "$OUTDIR/${sample}.pairs.gz") | \
samtools view -@ 8 -Sb - > "$OUTDIR/${sample}.bam"

echo "[5/6] Indexing BAM"
samtools index "$OUTDIR/${sample}.bam"

# Optional: Pairix indexing
# echo "[6/6] Indexing .pairs.gz"
# pairix "$OUTDIR/${sample}.pairs.gz"

rm "$OUTDIR/${sample}.dedup.pairsam"

echo " ^|^e All done for $sample"

