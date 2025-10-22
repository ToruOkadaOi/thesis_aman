#(hicexplorer_env) [papantonis1@gwdu101 cooltools_dots]$ wc -l both.txt only10.txt only5.txt
#  3346 both.txt
#  3382 only10.txt
#   832 only5.txt
#  7560 total

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

only10 = 3382
only5 = 832
both = 3346

venn2(subsets=(only10, only5, both), set_labels=("RBP1_10kb", "RBP1_5kb"))
plt.title("Overlap of RBP1 Loops in Merged Set (FDR < 0.01)")
plt.tight_layout()
plt.savefig("rbp1_loop_overlap_venn.png", dpi=300)
