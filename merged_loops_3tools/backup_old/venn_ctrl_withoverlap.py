import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import numpy as np

# ------------------------------
# Configuration
# ------------------------------
max_dist = 10000  # bp allowed between anchor midpoints
peakachu_file = "peakachu_merged_loops_ctrl.bedpe"
mustache_file = "mustache_merged_loops_ctrl.bedpe"
cooldots_file = "cooldots_merged_loops_ctrl.bedpe"

# ------------------------------
# Load and process loops
# ------------------------------
def load_loops_df(filepath):
    df = pd.read_csv(filepath, sep='\t', header=None, usecols=range(6))
    df.columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    df['mid1'] = ((df['start1'] + df['end1']) // 2).astype(int)
    df['mid2'] = ((df['start2'] + df['end2']) // 2).astype(int)
    # Normalize anchor order
    df[['a_chr', 'a_mid', 'b_chr', 'b_mid']] = df.apply(
        lambda row: (
            row['chr1'], row['mid1'], row['chr2'], row['mid2']
        ) if (row['chr1'], row['mid1']) <= (row['chr2'], row['mid2']) else (
            row['chr2'], row['mid2'], row['chr1'], row['mid1']
        ),
        axis=1, result_type='expand'
    )
    df['loop_id'] = df['a_chr'] + ':' + df['a_mid'].astype(str) + '-' + df['b_chr'] + ':' + df['b_mid'].astype(str)
    return df[['a_chr', 'a_mid', 'b_chr', 'b_mid', 'loop_id']]

# ------------------------------
# Extract the unique loops per tool as bedpe form
# ------------------------------

def extract_original_bedpe(df, id_set, out_file):
    df_filtered = df[df['loop_id'].isin(id_set)].copy()
    # Reconstruct original BEDPE columns
    df_filtered['start1'] = df_filtered['a_mid'] - 2500
    df_filtered['end1']   = df_filtered['a_mid'] + 2500
    df_filtered['start2'] = df_filtered['b_mid'] - 2500
    df_filtered['end2']   = df_filtered['b_mid'] + 2500
    df_filtered = df_filtered[['a_chr', 'start1', 'end1', 'b_chr', 'start2', 'end2']]
    df_filtered.to_csv(out_file, sep='\t', header=False, index=False)

# ------------------------------
# Fast fuzzy overlap via pandas
# ------------------------------
def fuzzy_overlap(df1, df2, max_dist):
    merged = df1.merge(df2, on=['a_chr', 'b_chr'], suffixes=('_1', '_2'))
    close = merged[
        (abs(merged['a_mid_1'] - merged['a_mid_2']) <= max_dist) &
        (abs(merged['b_mid_1'] - merged['b_mid_2']) <= max_dist)
    ]
    return set(close['loop_id_1'])

# ------------------------------
# Load loop sets
# ------------------------------
df_peakachu = load_loops_df(peakachu_file)
df_mustache = load_loops_df(mustache_file)
df_cooldots = load_loops_df(cooldots_file)

peakachu_ids = set(df_peakachu['loop_id'])
mustache_ids = set(df_mustache['loop_id'])
cooldots_ids = set(df_cooldots['loop_id'])

# ------------------------------
# Fuzzy overlaps
# ------------------------------
p_m = fuzzy_overlap(df_peakachu, df_mustache, max_dist)
p_c = fuzzy_overlap(df_peakachu, df_cooldots, max_dist)
m_c = fuzzy_overlap(df_mustache, df_cooldots, max_dist)

p_in_all = p_m & p_c

# Unique and pairwise
p_and_m_only = p_m - p_in_all
p_and_c_only = p_c - p_in_all
m_and_c_only = m_c - fuzzy_overlap(df_mustache[df_mustache['loop_id'].isin(m_c)], df_peakachu, max_dist)

only_p = peakachu_ids - (p_m | p_c)
only_m = mustache_ids - fuzzy_overlap(df_mustache, pd.concat([df_peakachu, df_cooldots]), max_dist)
only_c = cooldots_ids - fuzzy_overlap(df_cooldots, pd.concat([df_peakachu, df_mustache]), max_dist)

extract_original_bedpe(df_peakachu, only_p, "unique_peakachu_ctrl.bedpe")
extract_original_bedpe(df_mustache, only_m, "unique_mustache_ctrl.bedpe")
extract_original_bedpe(df_cooldots, only_c, "unique_cooldots_ctrl.bedpe")

# ------------------------------
# Plot Venn
# ------------------------------
venn3(
    subsets=(
        len(only_p),
        len(only_m),
        len(p_and_m_only),
        len(only_c),
        len(p_and_c_only),
        len(m_and_c_only),
        len(p_in_all)
    ),
    set_labels=('Peakachu', 'Mustache', 'Cooldots')
)

# Add total loop counts as annotations outside circles
plt.text(-0.9, 0.8, f"Peakachu total: {len(peakachu_ids)}", fontsize=10, color='black')
plt.text(0.9, 0.8, f"Mustache total: {len(mustache_ids)}", fontsize=10, color='black')
plt.text(0.0, -0.9, f"Cooltools-dots total: {len(cooldots_ids)}", fontsize=10, color='black', ha='center')

plt.title(f"Control loop overlap (±{max_dist // 1000}kb Midpoint)")
plt.savefig("mergedloops_ctrl.png", dpi=300, bbox_inches='tight')
plt.show()