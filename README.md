GENESIS: GENetic Effect-Size distribution Inference from Summary-level data
====
  
## Overview
To goal of GENESIS is to analyze summary-level GWAS statistics and external linkage disequilibrium information to estimate common variants effect-size distributions, characterized by the proportion of underlying susceptibility SNPs and a flexible normal-mixture model for their effects. This package allows model flexibility by considering a 2- or 3-component model, which means modeling the effect-size distributions for susceptibility SNPs by a single normal, or a mixture of two normal distributions. This package also allows users to make predictions of future GWAS with a bigger sample size.  

## GENESIS Installation

GENESIS software can be installed via Github. To install the latest version of GENESIS package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("GENESIS","yandorazhang")
```

## Citation

Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits and implications for the future." bioRxiv (2017): 175406.


## Contact the Author
Author: Yan Zhang, Guanghao Qi, Ju-Hyun Park, Nilanjan Chatterjee

Maintainer: Yan Zhang (yandorazhang@gmail.com)