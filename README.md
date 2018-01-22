GENESIS: GENetic Effect-Size distribution Inference from Summary-level data
====
  
## Overview

The goal of GENESIS is to analyze summary-level GWAS statistics and external linkage disequilibrium information to estimate common variants effect-size distributions, characterized by the proportion of underlying susceptibility SNPs and a flexible normal-mixture model for their effects. This package allows flexibility by considering a 2- or 3-component model, which respectively incorporate a single normal, or a mixture of two normal distributions, for specifying the effects of non-null SNPs. This package also allows users to make predictions regarding yield of future GWAS with larger sample sizes.

## GENESIS Installation

GENESIS software can be installed via Github. To install the latest version of GENESIS package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("yandorazhang/GENESIS")
```

## Citation

Please cite the following paper when you use GENESIS:

Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.


## Contact the Author
Author: Yan Zhang, Guanghao Qi, Ju-Hyun Park, Nilanjan Chatterjee

Maintainer: Yan Zhang (yandorazhang@gmail.com)
