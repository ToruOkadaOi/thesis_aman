#!/bin/bash
#SBATCH -p medium
#SBATCH --job-name=matrix_gen_clean
#SBATCH --output=logs/matrix_%A_%a.out
#SBATCH --error=logs/matrix_%A_%a.err
#SBATCH --array=0-23
#SBATCH --time=46:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1

# ----------------------------
# Set working directory
cd /usr/users/papantonis1/aman/microc_project

# Chromosome list (array task ID maps to this)
chromList=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
chrom=${chromList[$SLURM_ARRAY_TASK_ID]}

# Parameters
res_kb=5
matrix_size=21

# Data input path
DPATH="/usr/users/papantonis1/aman/microc_project/loop_calling_premade_hic/RBP1"

# Output path
save_path="/usr/users/papantonis1/aman/microc_project/loop_calling_premade_hic/chr_all_sample"
mkdir -p "$save_path"

# Output file names
tmp_file="${save_path}/chr${chrom}_matrixsize${matrix_size}_tmp.npy"
out_file="${save_path}/chr${chrom}_matrixsize${matrix_size}.npy"

# Absolute paths to scripts
CHR_SCRIPT="/home/uni08/papantonis1/aman/CGLoop/DataProcessing/chr_all_sample.py"
CTRL_SCRIPT="/home/uni08/papantonis1/aman/CGLoop/DataProcessing/control_contact.py"

# Run job
echo "[$(date)] Starting chr${chrom} at ${res_kb}kb"

python "$CHR_SCRIPT" \
  "${DPATH}/KR_matrix_${res_kb}kb.chr${chrom}" \
  "$tmp_file" \
  $matrix_size $res_kb

if [[ -f "$tmp_file" ]]; then
  python "$CTRL_SCRIPT" "$tmp_file" "$out_file"
  echo "[$(date)] Finished chr${chrom}"
else
  echo "[$(date)] ERROR: Failed to generate $tmp_file for chr${chrom}" >&2
fi

