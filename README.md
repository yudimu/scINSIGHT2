# scINSIGHT2
Harmonizing heterogeneous single-cell gene expression data with individual-level covariate information

## Introduction
We introduce scINSIGHT2, a new integration model designed to harmonize gene expression data from multiple single-cell samples by incorporating both discrete and continuous individual-level covariates. Utilizing a generalized linear latent variable model, scINSIGHT2 adjusts for covariate-associated gene expression changes prior to estimating cell embeddings within a unified low-dimensional space of inferred metagenes. 

## Installation
`scINSIGHT2` is currently available on GitHub. To download and load package please follow:

``` r
devtools::install_github("https://github.com/yudimu/scINSIGHT2")
library(scINSIGHT2)
```

The installation has been tested under macOS 15.4.1 with R version 4.2.1 (64-bit x86_64-apple-darwin17.0 platform). and Rocky Linux 8.8 (Green Obsidian) with R version 4.4.2 (64-bit x86_64-pc-linux-gnu platform). 

## Tutorial
Please see the Wiki page for a detailed [tutorial](https://github.com/yudimu/scINSIGHT2/wiki/scINSIGHT2-vignette) on how to use this package.

## Issues and communications

If you have any issues using this package, please post them
[here](https://github.com/yudimu/scINSIGHT2/issues). Any suggestions and
comments are welcome! For suggestions and comments, please contact
Yudi Mu (<ymu015@ucr.edu>) or Wei Vivian Li (<weil@ucr.edu>).
