<!-- README.md is generated from README.Rmd. Please edit that file -->
Quick start
-----------

Welcome to the `propr` GitHub page!

The bioinformatic evaluation of gene co-expression often begins with correlation-based analyses. However, this approach lacks statistical validity when applied to relative count data. This includes, for example, those biological data produced by microarray assays or high-throughput RNA-sequencing. This package provides a set of functions for measuring dependence between relative features using compositional data analysis. Specifically, this package implements two measures of *proportionality*, φ and ρ, introduced in Lovell 2015 and expounded in Erb 2016. You can get started with `propr` by installing the most up-to-date version of this package directly from GitHub.

``` r
library(devtools)
devtools::install_github("tpq/propr")
library(propr)
```

The principal functions in `propr` include: (1) `phit`, for the calculation of φ, and (2) `perb`, for the calculation of ρ. In the example below, we calculate these proportionality metrics for a simulated dataset, printing the results as a proportionality matrix and then as a pairwise summary. We refer you to the official package vignette for a comprehensive discussion of compositional data, proportionality, and everything this package has to offer.

``` r
set.seed(12345)
N <- 10
data.absolute <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
                            c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
data.relative <- data.absolute / colSums(data.absolute)
```

### Calculate φ

``` r
phi <- phit(data.relative)
```

    ## Calculating all phi for actual counts...

``` r
phi@matrix
```

    ##             a           b          c          d
    ## a 0.000000000 0.001894476 3.95056338 4.02312199
    ## b 0.001894476 0.000000000 3.97849497 4.05353543
    ## c 3.950563382 3.978494970 0.00000000 0.01119647
    ## d 4.023121991 4.053535432 0.01119647 0.00000000

``` r
phi@pairs
```

    ##   feature1 feature2        prop
    ## 1        a        b 0.001894476
    ## 2        c        d 0.011196465
    ## 3        a        c 3.950563382
    ## 4        b        c 3.978494970
    ## 5        a        d 4.023121991
    ## 6        b        d 4.053535432

### Calculate ρ

``` r
rho <- perb(data.relative)
```

    ## Calculating all rho for actual counts...

``` r
rho@matrix
```

    ##            a          b          c          d
    ## a  1.0000000  0.9990459 -0.9985539 -0.9982335
    ## b  0.9990459  1.0000000 -0.9981875 -0.9985699
    ## c -0.9985539 -0.9981875  1.0000000  0.9945048
    ## d -0.9982335 -0.9985699  0.9945048  1.0000000

``` r
rho@pairs
```

    ##   feature1 feature2       prop
    ## 1        a        b  0.9990459
    ## 2        b        d -0.9985699
    ## 3        a        c -0.9985539
    ## 4        a        d -0.9982335
    ## 5        b        c -0.9981875
    ## 6        c        d  0.9945048

### References

1.  Erb, I. & Notredame, C. 2016. How should we measure proportionality on relative gene expression data? Theory Biosci.

2.  Lovell, D. et al. 2015. Proportionality: A Valid Alternative to Correlation for Relative Data. PLoS Comput Biol 11.
