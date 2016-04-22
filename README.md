<!-- README.md is generated from README.Rmd. Please edit that file -->
Quick start
-----------

Welcome to the `propr` GitHub page!

The bioinformatic evaluation of gene co-expression often begins with correlation based analyses. However, this approach lacks statistical validity when applied to relative count data. This might include, for example, those biological data produced by microarray assays or high-throughput RNA-sequencing. This package provides a set of functions for evaluating co-expression between relative features using compositional data analysis. Specifically, this package implements two measures of *proportionality*, \(\phi\) and \(\rho\), introduced in Lovell 2015 and expounded in Erb 2016. You can get started with `propr` by installing the most up-to-date version of this package directly from GitHub.

``` r
library(devtools)
devtools::install_github("tpq/propr")
```

    ## Skipping install for github remote, the SHA1 (b1658bc8) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(propr)
```

The principal functions in `propr` include: (1) `phit`, for the calculation of \(\phi\), and (2) `perb`, for the calculation of \(\rho\). In the example below, we calculate these proportionality scores for simulated relative data, printing the results as a proportionality matrix and also pairwise. We refer you to the official vignette for a comprehensive discussion of compositional data, proportionality, and everything this package has to offer.

``` r
set.seed(12345)
N <- 10
X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
                c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
data.absolute <- t(X) # set features as rows...
data.relative <- data.absolute / colSums(data.absolute)
```

``` r
phi <- phit(data.relative)
```

    ## Calculating all phi for actual counts...

``` r
phi@matrix
```

    ##             a           b         c         d
    ## a 0.000000000 0.007348132 4.0301680 4.0941984
    ## b 0.007348132 0.000000000 3.9051763 3.9754729
    ## c 4.030168000 3.905176324 0.0000000 0.0192352
    ## d 4.094198403 3.975472865 0.0192352 0.0000000

``` r
phi@pairs
```

    ##   feature1 feature2        prop
    ## 1        a        b 0.007348132
    ## 2        c        d 0.019235204
    ## 3        b        c 3.905176324
    ## 4        b        d 3.975472865
    ## 5        a        c 4.030168000
    ## 6        a        d 4.094198403

``` r
rho <- perb(data.relative)
```

    ## Calculating all rho for actual counts...

``` r
rho@matrix
```

    ##            a          b          c          d
    ## a  1.0000000  0.9964378 -0.9979883 -0.9954194
    ## b  0.9964378  1.0000000 -0.9954814 -0.9980810
    ## c -0.9979883 -0.9954814  1.0000000  0.9905436
    ## d -0.9954194 -0.9980810  0.9905436  1.0000000

``` r
rho@pairs
```

    ##   feature1 feature2       prop
    ## 1        b        d -0.9980810
    ## 2        a        c -0.9979883
    ## 3        a        b  0.9964378
    ## 4        b        c -0.9954814
    ## 5        a        d -0.9954194
    ## 6        c        d  0.9905436

References
----------

\(^1\) Erb, Ionas, and Cedric Notredame. “How Should We Measure Proportionality on Relative Gene Expression Data?” Theory in Biosciences = Theorie in Den Biowissenschaften, January 13, 2016. <doi:10.1007/s12064-015-0220-8>.

\(^2\) Lovell, David, Vera Pawlowsky-Glahn, Juan José Egozcue, Samuel Marguerat, and Jürg Bähler. “Proportionality: A Valid Alternative to Correlation for Relative Data.” PLoS Computational Biology 11, no. 3 (March 16, 2015). <doi:10.1371/journal.pcbi.1004075>.
