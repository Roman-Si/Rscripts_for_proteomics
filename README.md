# Rsripts_for_proteomics

## About
This repository contains R functions to parse output files from proteomics pipelines and perform downstream exploratory and differential expression analysis.

Functions to do the following:
1. Parse `mzTab` output from **quantms**
2. Parse `msstats_input` files from **quantms**
3. Parse DIA-NN `parquet` output
4. Format dataframes for use with **[prolfqua](https://github.com/fgcz/prolfqua)**
5. Perform quantitative and differential abundance analysis with **[prolfqua](https://github.com/fgcz/prolfqua)**


Can add Rmd files or links to other projects later to show example usage.

---

## Installation

Install the required Bioconductor dependencies before installing the package via `devtools`.

```r
# Install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "MSnbase", "ProtGenerics"))

# Install devtools and this package
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Roman-Si/Rscripts_for_proteomics")
```

---

## Usage
Example analyses using functions from this package:

* Proteomic analysis of *Aspergillus parasiticus* MM36 secretomes grown on long-chain alkanes and hexadecane:
  [MS.html](https://htmlpreview.github.io/?https://github.com/Roman-Si/Asp-parasiticus_alkane_multiomics/blob/main/MS.html)

* Proteomic analysis of *Pleurotus pulmonarius* LGAM 28684 secretomes grown on xylose, corn stover, and beechwood:
  [prolfqua.html](https://htmlpreview.github.io/?https://github.com/Roman-Si/Pleurotus_proteomics/blob/main/docs/prolfqua.html)

---


## Template projects

### Proteomics R packages used:
* [prolfqua](https://github.com/fgcz/prolfqua)
* [MsnBase](https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html)

### Upstream pipelines
* [quantms](https://github.com/bigbio/quantms) for DDA
* [DIA-NN](https://github.com/vdemichev/DiaNN) for DIA


### Packaging inspiration
* [Building reproducible analytical pipelines with R](https://raps-with-r.dev/)
* [Reproducible and Trustworthy Workflows for Data Science](https://ubc-dsci.github.io/reproducible-and-trustworthy-workflows-for-data-science/)
