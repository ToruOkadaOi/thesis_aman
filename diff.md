# Differential Analysis — How Stable is Chromatin Architecture in the Absence of PRC2?

To assess how chromatin architecture responds to the loss of PRC2, I quantified changes in chromatin loops by merging loop sets from control and eed mutant DLD1 cells. The loops were classified as shared, gained, or lost using a 5 kb fuzzy overlap threshold (bedtools pairtopair, type=both).

Out of 32,277 total loops, 7,454 were shared between conditions, 13,241 were gained specifically in the eed mutant, and 11,582 were lost relative to control. This indicated substantial remodeling of the loop landscape (around 54% altered), though a core subset of loops remained preserved.

To quantify loop strength changes, I then extracted contact frequencies from the Hi-C contact matrices at 10 kb resolution and calculated log₂ fold changes (log₂FC = log₂[mutant signal / control signal]) for each shared loop. Then loops with log₂FC > 1 as upregulated (stronger in the mutant) and those with log₂FC < -1 as downregulated.

Among the shared loops, 52 were identified to be upregulated and 97 were downregulated in the eed mutant. Notably, 7 of the downregulated shared loops overlapped exactly with H3K27me3-enriched loops (n = 1183).

To specifically assess the impact of PRC2 on H3K27me3-associated loops, I examined how many loops had H3K27me3 peaks at both or either anchor. In control cells, 184 loops had H3K27me3 signal at both anchors, compared to only 87 in the eed mutant, indicating a 52.7% reduction. When considering loops with H3K27me3 at either anchor, I observed 1,183 in control and 643 in the mutant—a 45.6% decrease.

These data suggest that loss of PRC2 substantially diminishes the frequency of H3K27me3-associated loops. This supports the idea that PRC2 contributes to the maintenance of specific long-range contacts at Polycomb-repressed loci, while global loop architecture remains relatively intact.

I performed aggregate loop pileups over all H3K27me3-associated loops in control (n = 1,183) and eed mutant (n = 643) cells. Despite the ~45% reduction in loop number, the pileup signal remained comparably strong in the mutant (max intensity ≈ 4.44 vs. 4.4 in control), suggesting that the subset of loops retained upon PRC2 loss are structurally robust. These results indicate that while many H3K27me3-associated loops are lost, those that persist maintain similar contact strength, implying selective preservation of stable long-range contacts in the absence of PRC2.

I further compared the size distribution of H3K27me3-associated loops between control and eed mutant cells. In control cells, the mean loop size was 411 kb, while in the mutant it was 399 kb. Median values were similarly close (305 kb in control vs. 290 kb in mutant). The overall distribution—spanning from ~35–2870 kb in control and ~40–2820 kb in mutant—was nearly unchanged.

This suggests that although the total number of H3K27me3-associated loops decreased substantially upon PRC2 loss, the size range and structure of the remaining loops were preserved. Therefore, PRC2 does not define loop length per se, but likely stabilizes a subset of regulatory contacts within that typical size range.

These results suggest that chromatin architecture in the absence of PRC2 is partially stable at a global level but may exhibit targeted weakening.

---

## TADs

To evaluate topological domain organization, I applied RobusTAD to identify high-quality TADs in both control and eed mutant cells. The original Hq files contained 18,725 (control) and 17,811 (eed) raw domain calls. After filtering for high-confidence domains (Hq == 1 and score > 0.5), I retained 6,468 TADs in control and 5,457 in the eed mutant.

This represents a ~15.6% reduction in robust TADs following PRC2 loss. The decrease supports the interpretation that PRC2 contributes to the maintenance of topological insulation at select regions. But just raw numbers alone won’t do justice, here differential analysis is necessary. This was done using the tool diffDomain by using the original Hq files.

This approach classified domain changes into six types: single, complex, split, loss, merge, and zoom, reflecting different reorganization modes (Hua et al., 2024).

While RobusTAD identified 18,725 and 17,811 domains in control and eed mutant cells respectively, diffDomain returned 14,636 entries after pairwise alignment. This reduction reflects filtering for directly comparable domains across conditions, excluding unmatched or ambiguous boundaries. Therefore, the diffDomain output represents the subset of domains confidently classifiable into structural reorganization types such as split, merge, or complex.

Out of all domains, the majority (n = 9,613) were classified as single, indicating stable one-to-one correspondence between conditions. However, a substantial number of TADs showed complex reorganization (n = 2,730) or underwent splitting (n = 874), loss (n = 854), or merging (n = 565).

Subtype classification of single domains further revealed that most exhibited changes in strength (n = 9,279), while only a small fraction showed zoom effects (n = 334)—i.e., local domain shrinkage or compaction.

These results indicate that while a core of TADs remains structurally intact, a significant fraction undergoes reorganization, weakening, or reconfiguration in the absence of PRC2. This suggests that PRC2 contributes to domain insulation fidelity and hierarchical structure, with its loss driving local, but not wholesale, architectural remodeling.

---

## Compartments

To examine compartment-level chromatin organization, I assigned A/B compartments by the sign of the first eigenvector (E1) at 100 kb resolution. Using compartment label switches, I identified 3,087 compartment transitions in control cells and 3,181 in eed mutant cells, suggesting a marginal increase in compartment segmentation following PRC2 loss.

Compartment counts were nearly balanced in both conditions (CTRL: 1,544 A and 1,543 B; RBP1: 1,590 A and 1,590 B), indicating that global A/B structure remains largely unchanged.

These results suggest that compartment-level architecture is stable in the absence of PRC2, with no major shifts in A/B balance, supporting the conclusion that PRC2 loss induces local, not global, architectural disruptions.

To evaluate large-scale compartmental organization, I computed the first principal component (E1) from Hi-C contact matrices at 100 kb resolution and compared control and eed mutant (RBP1) conditions. The genome-wide Pearson correlation of E1 values was r = 0.987, indicating near-complete preservation of A/B compartment structure following PRC2 loss. This supports the conclusion that global compartment identity remains stable, and that PRC2-dependent architectural changes occur at finer scales (e.g., specific loops and TADs).

I am considering a differential compartment analysis tool to implement ASAP.

To determine how compartment identity influences loop stability upon PRC2 loss, I annotated chromatin loops based on the A/B compartment status of their anchors, classifying them into types such as B_to_A_B_to_A, B_to_A_stable, and stable_B_to_A. These labels reflect whether loop anchors transitioned from inactive (B) to active (A) compartments or remained stable across conditions.

Among the mutant-specific loops, B_to_A_B_to_A loops were the most frequent (n = 45), followed by stable_B_to_A (n = 12) and B_to_A_stable (n = 11). A similar trend was observed among shared loops, where B_to_A_B_to_A loops (n = 24) remained the most common, followed by B_to_A_stable (n = 13) and stable_B_to_A (n = 10). This indicates that even in the absence of PRC2, loops spanning compartment-switching regions are not only formed de novo but also retained, suggesting that compartment reprogramming does not preclude stable long-range interactions.

The presence of shared loops with mixed or dynamic compartment identities implies that PRC2 is not strictly required for the maintenance of loops that bridge transcriptionally distinct regions. These data further support a model in which PRC2 selectively stabilizes a subset of chromatin contacts, particularly those associated with Polycomb repression, while leaving structurally persistent or compartmentally plastic interactions largely intact.

While A/B compartments are large-scale chromatin domains (typically 100 kb–1 Mb), loop anchors were annotated by the compartment identity of the 100 kb bin they reside in. This approach captures local compartmental context, rather than inferring that a loop necessarily spans entire compartments. The enrichment of loop types such as B_to_A_B_to_A and B_to_A_stable among shared and mutant-specific loops suggests that local transitions in chromatin environment do not prevent the formation or maintenance of long-range contacts, particularly in the absence of PRC2.

Would be interesting to see the interpretations from a published tool.

---

## Stripes

Differential stripe analysis
