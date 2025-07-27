<!-- Removed top-level header to avoid duplicate title -->
*Rendered with Cayman* *| Work in progress*

## Project Title

**Genome-wide profiling of H3K27me3 and its functional integration with 3D genome architecture in human cells**

### Location

[**Papantonis Lab**](https://papantonislab.eu), *Universitätsmedizin Göttingen*

## Background

To be updated

#### Data

- **Micro-C, CUT&Tag, RNA-Seq**: In-house generated datasets  
- **ATAC-Seq**: Public dataset from Lyu, Huijue, et al.  
  *“GATA6 is a novel regulator of gene expression and 3D genome in colorectal cancer.”*  
  *Cancer Research 84.6_Supplement (2024): 4401–4401*  
- **ChIP-Seq**: Retrieved from [ChIP-Atlas](https://chip-atlas.org) using query: `DLD-1`

## Chapters

> - Chapter.1 - [Chromatin Architecture in Wild-Type and EED-Mutant Conditions](loops.md)
> - Chapter.2 - [Differential analysis - How stable is chromatin architecture in the absence of PRC2?](diff.md)
> - Chapter.3 - [Integrating Chromatin Architecture with Transcriptional Dynamics](deg.md)
> - Chapter.4 - [Methodological Framework: Tools and Optimization](tools.md)
> -  Chapter.5 - [Exploratory Analyses](aux.md)


### Micro-C and Loop Analysis

| Category                  | Tool(s) / Documentation |
|--------------------------|-------------------------|
| QC & .pairs generation   | [Documentation](https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html) |
| Loop calling             | [Mustache](https://github.com/ay-lab/mustache), [Cooltools-dots](https://cooltools.readthedocs.io/en/latest/notebooks/dots.html), [Peakachu](https://github.com/tariks/peakachu) |
| TAD calling              | [RobusTAD](https://github.com/zhyanlin/RobusTAD), [Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead) |
| A/B compartment analysis | [Cooltools eigs-cis](https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html), [Calder2](https://github.com/CSOgroup/CALDER2) |
| Stripes                  | [Stripenn](https://github.com/ysora/stripenn), [Stripecaller](https://github.com/XiaoTaoWang/StripeCaller) |
| Differential analysis    | [diffDomain](https://github.com/Tian-Dechao/diffDomain) |

### RNA-Seq Analysis

| Category                  | Tool(s) |
|--------------------------|---------|
| Gene body analysis       | DESeq2 |
| Intron & Exon            | [iRNASeq](https://academic.oup.com/nar/article/43/6/e40/2453418?login=true) |
| GO & Genomic coordinates | [clusterProfiler](https://bioconductor.org/packages/devel/bioc/html/clusterProfiler.html), [TxDb](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html), [MyGene.py](https://github.com/biothings/mygene.py) |


## ML Modeling

| Task          | Model(s)                            | Purpose                                        |
|---------------|--------------------------------------|------------------------------------------------|
| Classification| Logistic Regression(baseline), XGBoost        | Classify loop status (gained, lost, shared)    |
| Regression     | XGBoost, LightGBM, CatBoost         | Predict gene expression (log₂FC) from features |


## Questions

In case of questions about the code or the study, please reach out [here](aman.nalakath@tum.de)

## References

## Acknowledgments

Many thanks to everyone at Papantonis lab. Also grateful for support from [Johannes lab](https://www.johanneslab.org/home)

### TO DOs

- HiGlass server 
- Nextflow for ML feature table generation
- [Ablation testing](https://medium.com/@rajilini/ablation-study-what-is-it-in-machine-learning-0a1d362b366d) in models

---
```
 Turns out PRC2 wasn't repressing genes…it was just trying to exploit some chromatin loop holes.
```
# ![Loophole High Five](./brooklyn-nine-nine-loophole.gif)
