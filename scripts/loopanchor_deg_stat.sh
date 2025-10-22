#!/bin/bash
# --------------------------------------------------
# DEG TSS vs loop anchors: proximity + distance tests
# Runs both UP and DOWN in one go
# --------------------------------------------------

# Paths
outdir=/usr/users/papantonis1/aman/rnaseq_data
anchors=$outdir/anchors.merged.bed

# Input TSS files
up_tss=$outdir/up_TSS_DEGs_gtfgencode.bed
bg_up_tss=$outdir/bg_up_TSS.bed
down_tss=$outdir/down_TSS_DEGs_gtfgencode.bed
bg_down_tss=$outdir/bg_down_TSS.bed

# --- Functions ---

run_proximity () {
    label=$1; tss=$2; bg=$3
    for w in 5000 10000; do
        bedtools window -a $tss -b $anchors -w $w \
          | cut -f1-3 | sort -u > $outdir/${label}_near_${w}.bed
        bedtools window -a $bg -b $anchors -w $w \
          | cut -f1-3 | sort -u > $outdir/bg_${label}_near_${w}.bed

        n_total=$(wc -l < $tss)
        n_near=$(wc -l < $outdir/${label}_near_${w}.bed)
        n_bg_total=$(wc -l < $bg)
        n_bg_near=$(wc -l < $outdir/bg_${label}_near_${w}.bed)

        python - "$n_total" "$n_near" "$n_bg_total" "$n_bg_near" "$w" "$label" <<'PY'
import sys
from scipy.stats import fisher_exact
n_total, n_near, n_bg_total, n_bg_near, w, label = sys.argv[1:]
n_total, n_near, n_bg_total, n_bg_near, w = map(int, [n_total, n_near, n_bg_total, n_bg_near, w])
table = [[n_near, n_total-n_near],
         [n_bg_near, n_bg_total-n_bg_near]]
odds, p = fisher_exact(table)
print(f"[Proximity] {label.upper()} | Window {w} bp | 2x2={table} | OR={odds:.2f} | p={p:.3e}")
PY
    done
}

run_distance () {
    label=$1; tss=$2; bg=$3
    bedtools closest -d -a $tss -b $anchors > $outdir/${label}_to_anchor.closest.tsv
    bedtools closest -d -a $bg  -b $anchors > $outdir/bg_${label}_to_anchor.closest.tsv

    python - "$outdir/${label}_to_anchor.closest.tsv" "$outdir/bg_${label}_to_anchor.closest.tsv" "$label" <<'PY'
import sys, pandas as pd, numpy as np
from scipy.stats import mannwhitneyu

f1, f2, label = sys.argv[1:]
df1 = pd.read_csv(f1, sep="\t", header=None)
df2 = pd.read_csv(f2, sep="\t", header=None)
d1 = df1.iloc[:,-1].astype(int)
d2 = df2.iloc[:,-1].astype(int)

print(f"\n[Distance] {label.upper()}")
print(f"Median distance ({label.upper()}): {np.median(d1):.0f} bp")
print(f"Median distance (BG): {np.median(d2):.0f} bp")
u,p = mannwhitneyu(d1, d2, alternative='two-sided')
print(f"Mann–Whitney p={p:.3e} | Δ-median={np.median(d2)-np.median(d1):.0f} bp")
PY
}

# --- Run for UP ---
run_proximity "up"   $up_tss   $bg_up_tss
run_distance  "up"   $up_tss   $bg_up_tss

# --- Run for DOWN ---
run_proximity "down" $down_tss $bg_down_tss
run_distance  "down" $down_tss $bg_down_tss

