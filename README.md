# Rsripts_for_proteomics

## About
Sripts for parsing the output of proteomics analyses. Started writing for the mzTab and msstats_input files from quantms (https://doi.org/10.1038/s41592-024-02343-1) but can expand.

Can add Rmd files or links to other projects later to show example usage.

## Installation

```
# if devtools not installed
# install.packages("devtools")

library(devtools)
devtools::install_github("Roman-Si/Rsripts_for_proteomics")
```

## Usage
An analysis using this package can be found in https://github.com/Roman-Si/Asp-parasiticus_alkane_multiomics 


## Template projects

Proteomics R packages used:
- [prolfqua](https://github.com/fgcz/prolfqua)
- [MsnBase](https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html)


Ideas for packaging taken from
- Building reproducible analytical pipelines with R, https://raps-with-r.dev/
- Reproducible and Trustworthy Workflows for Data Science, https://ubc-dsci.github.io/reproducible-and-trustworthy-workflows-for-data-science/