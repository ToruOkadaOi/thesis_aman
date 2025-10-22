#!/bin/bash
#SBATCH --job-name=xgboost_mc
#SBATCH --output=logs/xgboost_%j.out
#SBATCH --error=logs/xgboost_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=medium

# === Activate environment ===
source ~/.bashrc
conda activate cool_notebook

# === Run script ===
python /usr/users/papantonis1/aman/microc_project/loop_calling_premade_hic/jupyter_notes/mlmodel_improved.py
echo Done
