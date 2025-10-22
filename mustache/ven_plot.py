import matplotlib
matplotlib.use('Agg')  # for HPCs without display

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def load_loops(filepath):
    df = pd.read_csv(filepath, sep='\t', header=None, usecols=[0,1,2,3,4,5])
    return set(df.apply(lambda x: "{}:{}-{}_{}:{}-{}".format(x[0], x[1], x[2], x[3], x[4], x[5]), axis=1))

# Load 3 filtered loop sets
s1 = load_loops("RBP1_10kb_fdr001_clean.bedpe")
s2 = load_loops("RBP1_5kb_fdr001_clean.bedpe")
s3 = load_loops("merged_loops_rbp1.bedpe")

# Plot
venn3([s1, s2, s3], set_labels=('RBP1_10kb', 'RBP1_5kb', 'merged'))
plt.title("Loop Overlap (FDR < 0.01)")
plt.tight_layout()
plt.savefig("venn_fdr001_rbp1.png", dpi=300)

