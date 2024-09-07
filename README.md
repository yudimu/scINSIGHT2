# scINSIGHT2
Integration of heterogeneous single-cell gene expression data including individual-level information

## Introduction
We proposed an integration method called scINSIGHT2, which takes an expression count matrix and individual-level covariates (such as age, gender, etc.) as inputs, and returns latent factors representing the integrated data for downstream analyses such as clustering. The novelty of our method lies in its ability to consider both continuous and discrete individual-level covariates.

## Installation
`scINSIGHT2` is currently available on GitHub. To download and load package please follow:

``` r
devtools::install_github("https://github.com/yudimu/scINSIGHT2")
library(scINSIGHT2)
```

## Issues and communications

If you have any issues using this package, please post them
[here](https://github.com/yudimu/scINSIGHT2/issues). Any suggestions and
comments are welcome! For suggestions and comments, please contact
Yudi Mu (<ymu015@ucr.edu>) or Wei Vivian Li (<weil@ucr.edu>).
