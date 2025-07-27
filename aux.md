## Differential Contact Maps Based on Gene Expression Changes

Is there a tool that can generate a **differential chromatin contact map** based on expression changes, given that a control condition and its corresponding contact matrix are available?

I'm exploring whether any existing tool can integrate gene expression changes with chromatin contact data to generate **condition-specific or differential contact maps**

---

## CLOD Reimplementation
todo - add the paper here

To identify **long-range chromatin contacts** exhibiting condition-specific enrichment, I reimplemented the **CLOD** (Chromosome Long-range colocalization identifier through Outlier Detection) method from *Zhang et al., 2024*.

The method operates as follows:

- Smooths the Hi-C contact matrix using a **uniform filter**
- Computes **genomic distances** between all contact pairs
- Groups contacts into bins based on separation distance
- Within each distance bin:
  - Computes the **third quartile (Q3)** and **interquartile range (IQR)**
  - Flags contacts as outliers if:
    ```
    Contact signal > Q3 + 2 × IQR
    ```

This outlier-based detection is designed to isolate **statistically enriched, high-contact interactions**.

I applied this pipeline at **10 kb resolution** on the **balanced Hi-C matrix of chromosome 14**. 

> This minimal yet functional reimplementation captures the statistical foundation of CLOD

---

## Evaluation of CLOD-Detected Contacts

While CLOD detects statistical outliers, the corresponding **Aggregate Peak Analysis (APA)** plot shows **low to modest central enrichment**.

- APA reflects the **mean behavior** of the contact set
- Signal **heterogeneity** across contacts(esp. long distance) may blunt the APA signal
- Variable interaction patterns can reduce sharpness despite real biological significance

---

## Additional Analyses

### Motif Enrichment  
(todo)

### GC Content  
(todo)

