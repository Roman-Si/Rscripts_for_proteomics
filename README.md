# Rsripts_for_proteomics

## About
Scripts for parsing the output of proteomics analyses.

Functions to do the following:
1. Parse the mzTab output from quantms
2. Parse the msstats_input file from quantms
3. Parse the diann parquet output
4. Prepare dataframes for EDA and differential analysis with prolfqua

Can add Rmd files or links to other projects later to show example usage.

## Installation

You can install with devtools, however devtools does not resolve the Bioconductor packages. They have to be installed beforehand.

```
# Install Bioconductor packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "MSnbase", "ProtGenerics"))

# Install devtools and this package
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
devtools::install_github("Roman-Si/Rsripts_for_proteomics")
```

## Usage
An analysis using this package can be found in https://github.com/Roman-Si/Asp-parasiticus_alkane_multiomics 


## Template projects

Proteomics R packages used:
- [prolfqua](https://github.com/fgcz/prolfqua)
- [MsnBase](https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html)

Pipelines used for analysis of LC-MS/MS data:
- [quantms](https://github.com/bigbio/quantms) for DDA
- [DIA-NN](https://github.com/vdemichev/DiaNN) for DIA


Ideas for packaging taken from
- Building reproducible analytical pipelines with R, https://raps-with-r.dev/
- Reproducible and Trustworthy Workflows for Data Science, https://ubc-dsci.github.io/reproducible-and-trustworthy-workflows-for-data-science/
