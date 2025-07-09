import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

def load_ids(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

# Load loop ID sets for CTRL
peak = load_ids("exact_ids/peakachu_merged_loops_ctrl_ids.txt")
must = load_ids("exact_ids/mustache_merged_loops_ctrl_ids.txt")
cool = load_ids("exact_ids/cooldots_merged_loops_ctrl_ids.txt")

# Plot Venn diagram
venn3([peak, must, cool], set_labels=["Peakachu", "Mustache", "Cooltools"])
plt.title("CTRL Loop Overlap (Exact Match)")
plt.tight_layout()
plt.savefig("venn_ctrl_exact.png", dpi=300)