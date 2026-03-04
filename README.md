# ***rnaseqtools***

**Downstream RNA-Seq Data Analysis and Visualization**

A collection of tools for RNA-seq data analysis, including feature count manipulation and summarization, gene set testing and data visualization.


---


## Description

`rnaseqtools` provides streamlined functions for:
- downloading and processing GENCODE annotations and MSigDB collections to produce enriched feature metadata.
- manipulating, transforming and summarizing feature counts (e.g., TPM rescaling, pseudobulk computation, PCA).
- standardizing differential gene expression results produced by the `DESeq2` and `MAST` packages.
- running gene set tests from the `fgsea` and `limma` packages on multiple gene collections.
- data visualization (e.g., plotting PCA, DGEA and GSEA results).


---


## Installation

You can install the development version of `rnaseqtools` from [GitHub](https://github.com/) with:

```r
remotes::install_github("benvallin/rnaseqtools")
```


---
