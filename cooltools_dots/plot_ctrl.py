#(hicexplorer_env) [papantonis1@gwdu101 cooltools_dots]$ wc -l both_ctrl.txt only10_ctrl.txt only5_ctrl.txt
#  3838 both_ctrl.txt
#  3671 only10_ctrl.txt
#  1045 only5_ctrl.txt
#  8554 total

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

only10 = 3671
only5 = 1045
both = 3838

venn2(subsets=(only10, only5, both), set_labels=("ctrl_10kb", "ctrl_5kb"))
plt.title("Overlap of ctrl Loops in Merged Set (FDR < 0.01)")
plt.tight_layout()
plt.savefig("ctrl_loop_overlap_venn.png", dpi=300)
