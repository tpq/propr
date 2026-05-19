<!-- README.md is generated from README.Rmd. Please edit that file -->

## Attention!!

If you have used `propr` previously, you will notice some changes. In
version 5.0.0, I have completed a major revision of the code base to
simplify the maintenance of `propr` going forward, including a
restructure of back-end and front-end functionality. From the user’s
perspective, you will notice a few changes. First, all supporting
visualization functions are gone. They were poorly implemented and not
backwards compatible. Instead, you can use the unified `getResults`
wrapper to pull data from `propr` and `propd` objects to pipe to
`ggplot2` for visualization. Second, many experimental functions have
been removed, with the remaining ones all sharing the prefix `run`. The
core routines called by the `propr` and `propd` functions remain
unchanged.

## Requirements

- GCC 13.3.0 (doesn't compile yet on GCC 14)

## Introduction

The `propr` package provides an interface for 4 distinct approaches to
compositional data analysis (CoDA): proportionality, differential
proportionality, logratio partial correlation with basis shrinkage, and
ratio analysis.

If you use this software, please cite our work. We don’t get paid to
make software, but your citations help us to negotiate support for
software maintenance and development.

``` r
citation("propr")
```

    ## 
    ## To cite propr in publications use:
    ## 
    ##   Jin S, Notredame C, Erb I (2023) Compositional covariance shrinkage
    ##   and regularised partial correlations. Statistics and Operations
    ##   Research Transactions 47(2). doi:10.57645/20.8080.02.8
    ## 
    ##   Quinn TP, Erb I, Gloor G, Notredame C, Richardson MF, Crowley TM
    ##   (2019) A field guide for the compositional analysis of any-omics
    ##   data. GigaScience 8(9). doi:10.1093/gigascience/giz107
    ## 
    ##   Quinn T, Erb I, Richardson MF, Crowley T (2018) Understanding
    ##   sequencing data as compositions: an outlook and review.
    ##   Bioinformatics 34(16): doi:10.1093/bioinformatics/bty175
    ## 
    ##   Erb I, Quinn T, Lovell D, Notredame C (2017) Differential
    ##   Proportionality - A Normalization-Free Approach To Differential Gene
    ##   Expression. Proceedings of CoDaWork 2017, The 7th Compositional Data
    ##   Analysis Workshop; available under bioRxiv 134536: doi:10.1101/134536
    ## 
    ##   Quinn T, Richardson MF, Lovell D, Crowley T (2017) propr: An
    ##   R-package for Identifying Proportionally Abundant Features Using
    ##   Compositional Data Analysis. Scientific Reports 7(16252):
    ##   doi:10.1038/s41598-017-16520-0
    ## 
    ##   Erb I, Notredame C (2016) How should we measure proportionality on
    ##   relative gene expression data? Theory in Biosciences 135(1):
    ##   doi:10.1007/s12064-015-0220-8
    ## 
    ##   Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S, Bahler J (2015)
    ##   Proportionality: A Valid Alternative to Correlation for Relative
    ##   Data. PLoS Computational Biology 11(3):
    ##   doi:10.1371/journal.pcbi.1004075
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.

OK, now let’s get started.

``` r
counts <- matrix(rpois(20*50, 100), 20, 50)
group <- sample(c("A", "B"), size = 20, replace = TRUE)
devtools::install_github("tpq/propr")
library(propr)
```

## Proportionality

There are a few proportionality statistics available. Select one with
the ‘metric’ argument.

``` r
pr <- propr(
        counts,  # rows as samples, like it should be
        metric = "rho",  # or "phi", "phs", "vlr"
        ivar = "clr",  # or can use another gene as reference, by giving the name or index
        alpha = NA,  # use to handle zeros
        p = 100  # used for permutation tests
      ) 
```

You can determine the “signficance” of proportionality using a built-in
permutation procedure. It estimates the false discovery rate (FDR) for
any cutoff. This method can take a while to run, but is parallelizable.

``` r
pr <- updateCutoffs(
        pr,
        number_of_cutoffs = 100,  # number of cutoffs to estimate FDR
        custom_cutoffs = NA,  # or specify custom cutoffs
        tails = 'right',  # consider only the positive values ('right') or both sides ('both')
        ncores = 4  # parallelize here
      ) 
```

Choose the largest cutoff with an acceptable FDR.

## Logratio partial correlation with basis shrinkage

There are many ways to calculate partial correlations, with or without
shrinkage. The recommended one for datasets with p\>\>n and influenced
by compositional bias is “pcor.bshrink”.

``` r
pr <- propr(
        counts,  # rows as samples, like it should be
        metric = "pcor.bshrink",  # partial correlation without shrinkage "pcor" is also available
        ivar = "clr",  # "clr" is recommended
        p = 100  # used for permutation tests
      ) 
```

You can also determine the “significance” of logratio partial
correlations with the built-in permutation approach.

``` r
pr <- updateCutoffs(
        pr,
        number_of_cutoffs = 100,  # number of cutoffs to estimate FDR
        custom_cutoffs = NA,  # or specify custom cutoffs
        tails = 'right',  # consider only the positive values ('right') or both sides ('both')
        ncores = 4  # parallelize here
      ) 
```

## Differential Proportionality

There are also a few differential proportionality statistics, but they
all get calculated at once.

``` r
pd <- propd(
        counts,
        group,  # a vector of 2 or more groups
        alpha = NA,  # whether to handle zeros
        p = 100,  # used for permutation tests
        weighted = TRUE  # whether to weight log-ratios
      )
```

You can switch between the “disjointed” and “emergent” statistics.

``` r
setDisjointed(pd)
```

``` r
setEmergent(pd)
```

You can again permute an FDR with the `updateCutoffs` method.
Alternatively, you can calculate an exact p-value for *θ* based on a
F-test. This is handled by the `updateF` method.

``` r
pd <- updateF(
        pd,
        moderated = FALSE,  # moderate stats with limma-voom
        ivar = "clr"  # used for moderation
      ) 
```

## Getters

Both functions return S4 objects. This package includes several helper
functions that work for both the `propr` and `propd` output.

``` r
?getMatrix # get results as a square matrix
?getResults # get propr or propd results in long-format
?getRatios # get samples by ratios matrix
```

Use `getResults` to pipe to `ggplot2` for visualization.

We also provide accesory functions to get the significant pairs.

``` r
?getSignificantResultsFDR
?getSignificantResultsFstat
?getAdjacencyFDR
?getAdjacencyFstat
?getCutoffFDR
?getCutoffFstat
```

Notice that for the getter functions to work properly on `propd`
objects, you have to set the target `theta` value active:

``` r
setActive(pd, "theta_d")
setActive(pd, "theta_e")
setActive(pd, "theta_mod")
```

## Ratio Methods

COMING SOON!!



## For Local Developement
For local developement install the pre-push script `./scripts/git/pre-push` to run the tests.