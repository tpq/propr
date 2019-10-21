<!-- README.md is generated from README.Rmd. Please edit that file -->
The `propr` package implements two analyses: proportionality and differential proportionality.

Proportionality is a compositional alternative to correlation. It is introduced by Lovell et al., expounded by Erb & Notredame, and implemented by Quinn et al. (2017). Differential proportionality is a compositional alternative to differential abundance. It is introduced by Erb et al. (2017), and discussed further by Quinn et al. (2019). Compositional data analysis for genomics is reviewed by Quinn et al. (2018).

If you use this software, please cite our work. We don't get paid to make software, but your citations allow us to negotiate time allocation for software maintence and development.

``` r
citation("propr")
```

    ## 
    ## To cite propr in publications use:
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
    ##   Proportionality - A Normalization-Free Approach To Differential
    ##   Gene Expression. Proceedings of CoDaWork 2017, The 7th
    ##   Compositional Data Analysis Workshop; available under bioRxiv
    ##   134536: doi:10.1101/134536
    ## 
    ##   Quinn T, Richardson MF, Lovell D, Crowley T (2017) propr: An
    ##   R-package for Identifying Proportionally Abundant Features Using
    ##   Compositional Data Analysis. Scientific Reports 7(16252):
    ##   doi:10.1038/s41598-017-16520-0
    ## 
    ##   Erb I, Notredame C (2016) How should we measure proportionality
    ##   on relative gene expression data? Theory in Biosciences 135(1):
    ##   doi:10.1007/s12064-015-0220-8
    ## 
    ##   Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S, Bahler J
    ##   (2015) Proportionality: A Valid Alternative to Correlation for
    ##   Relative Data. PLoS Computational Biology 11(3):
    ##   doi:10.1371/journal.pcbi.1004075
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.

OK, now let's get started.

``` r
counts <- matrix(rpois(20*50, 100), 20, 50)
group <- sample(c("A", "B"), size = 20, replace = TRUE)
devtools::install_github("tpq/propr")
library(propr)
```

Proportionality
---------------

There are a few proportionality statistics available. Select one with the 'metric' argument.

``` r
pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
```

You can determine the "signficance" of proportionality using a built-in permutation procedure. It tells estimates the false discovery rate (FDR) for any cutoff. This method can take a while to run, but is parallelizable.

``` r
updateCutoffs(pr,
              cutoff = seq(0, 1, .05), # cutoffs at which to estimate FDR
              ncores = 1) # parallelize here
```

Choose the largest cutoff with an acceptable FDR.

Differential Proportionality
----------------------------

There are also a few differential proportionality statistics, but they all get calculated at once.

``` r
pd <- propd(counts,
            group, # a vector of 2 or more groups
            alpha = NA, # whether to handle zeros
            weighted = TRUE, # whether to weigh log-ratios
            p = 100) # used by updateCutoffs
```

You can switch between the "disjointed" and "emergent" statistics.

``` r
setDisjointed(pd)
```

``` r
setEmergent(pd)
```

You can again permute an FDR with the `updateCutoffs` method. Alternatively, you can calculate an exact p-value for *θ* based on a F-test. This is handled by the `updateF` method.

``` r
pd <- updateF(pd,
              moderated = FALSE, # moderate stats with limma-voom
              ivar = "clr") # used for moderation
```

Getters
-------

Both functions return S4 objects. This package includes several helper functions that work for both the `propr` and `propd` output. Most of the time, you would want to use `getResults`. This method only selects pairs beyond a certain size, chosen by the 'cutoff' argument.

``` r
?getResults # get results in long-format
?getMatrix # get results as a square matrix
?getAdj # get an adjacency matrix
```

The vignettes describe some custom visualization methods.
