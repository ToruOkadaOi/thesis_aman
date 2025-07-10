<!-- Removed top-level header to avoid duplicate title -->
*Rendered with Cayman*

## Project Title and Overview
Genome-wide analysis of H3K27me3 distribution and its functional correlation with 3D genome
architecture in human cells
### *Translational Epigenomics in 3D Genome Architecture*
### Work done at [Papantonis lab](https://papantonislab.eu), Universitätsmedizin Göttingen

Integration of Micro-C & RNA-Seq to understand loop dynamics under PRC2 inhibition.

## Background

To be updated

#### Data

## Chapters

- Chapter.1 - [Architectural changes](loops.md)

## Micro-C and Loop Analysis
- QC & .pairs file generation - [documentation](https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html)
- Loop calling via [Mustache](https://github.com/ay-lab/mustache), [Cooltools-dots](https://cooltools.readthedocs.io/en/latest/notebooks/dots.html), [Peakachu](https://github.com/tariks/peakachu)
- TAD calling - [RobusTAD](https://github.com/zhyanlin/RobusTAD) & [Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead)
- A/B compartment analysis - [Cooltools eigs-cis](https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html), [Calder2](https://github.com/CSOgroup/CALDER2)
- Stripes - [Stripenn](https://github.com/ysora/stripenn), [Stripecaller](https://github.com/XiaoTaoWang/StripeCaller)
- Differential analysis - [diffDomain](https://github.com/Tian-Dechao/diffDomain)

## RNA-Seq Analysis
Gene body analysis - DESeq2\
Intron & Exon - [iRNASeq](https://academic.oup.com/nar/article/43/6/e40/2453418?login=true)\
GO & Genomic coordinates - clusterProfiler, GSEA.py & TxDb

## ML Modeling
Logistic Regression & Gradient Boost Model

## Contact

aman.nalakath@tum.de or amanshamil@protonmail.com

### to-dos

#### HiGlass server for visualizing loops - 
to be added...

# ![Loophole High Five](./brooklyn-nine-nine-loophole.gif)
