<script>
  const link = document.createElement('link');
  link.rel = 'shortcut icon';
  link.type = 'image/x-icon';
  link.href = '/favicon.ico';
  document.head.appendChild(link);
</script>

<!-- Removed top-level header to avoid duplicate title -->
*Rendered with Cayman* *| Work in progress*

## Project Title and Overview

**Genome-wide analysis of H3K27me3 distribution and its functional correlation with 3D genome architecture in human cells**

### *Translational Epigenomics in 3D Genome Architecture*

**Work done at** [**Papantonis Lab**](https://papantonislab.eu), *Universitätsmedizin Göttingen*

Integration of Micro-C & RNA-Seq to understand loop dynamics under PRC2 inhibition.

## Background

To be updated

#### Data

## Chapters

- Chapter.1 - [Architectural changes](loops.md)
- Chapter.2 - [Differential analysis - How stable is chromatin architecture in the absence of PRC2?](diff.md)
- Chapter.3 - [Integrating Chromatin Architecture with Transcriptional Dynamics](deg.md)
- Chapter.4 - [Methodological Framework: Tool Performance and Optimization](tools.md)
- Chapter.5 - [Exploratory Analyses and Ancillary Discoveries](aux.md)

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
| GO & Genomic coordinates | clusterProfiler, GSEA.py, TxDb |


## ML Modeling

| Task          | Model(s)                            | Purpose                                        |
|---------------|--------------------------------------|------------------------------------------------|
| Classification| Logistic Regression(baseline), XGBoost        | Classify loop status (gained, lost, shared)    |
| Regression     | XGBoost, LightGBM, CatBoost         | Predict gene expression (log₂FC) from features |


## Questions

In case of questions about the code or the study, please reach out [here](aman.nalakath@tum.de)

## References

### TO DOs

HiGlass server for visualizing loops - 
to be added...

# ![Loophole High Five](./brooklyn-nine-nine-loophole.gif)
