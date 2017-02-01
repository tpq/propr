## propr 2.1.2
---------------------
* Modified visualization tools
    * New `pca` function masks `mds` function
* Include inst/CITATION file

## propr 2.1.1
---------------------
* Modified `propr` Class
    * Merged `propr-class` and `propr` documentation
* Modified `phit`, `perb` functions
    * Merged `phit` and `perb` documentation
    * New `phis` function returns `(1 - rho) / (1 + rho)`
    * NAs in count matrix now throw error
    * 0s now replaced with 1s
* Modified visualization tools
    * Merged documentation

## propr 2.1.0
---------------------
* Modified `[` method
    * Now joins newly indexed pairs with any existing index
* New `cytescape` function
    * Uses `@pairs` slot to build an interaction network

## propr 2.0.9
---------------------
* Modified visualization tools
    * Courtesy `prompt` argument extended to `smear` and `dendrogram`
    * Improved error handling and documentation
* Modified `abstract` function
    * New `dt` argument indexes significant results in `@pairs`
* Modified `simplify` function
    * Now builds index of lower left triangle of matrix
* New `adjacent` function
    * Uses `@pairs` slot to build an adjacency matrix

## propr 2.0.8
---------------------
* Modified visualization tools
    * `bucket` now depends on `slate` function
* Modified backend code
    * New `coordToIndex` performs inverse of `indexToCoord`
* Modified `prop2prob` function
    * Return p-values as a sorted `data.table`
    * Now lets user select `method` for p-value adjustment
    * New `prompt` argument turns off big data prompt
    * Fix pass by reference bug in `linRcpp`
* New `abstract` function
    * Combines two `propr` objects into one

## propr 2.0.7
---------------------
* New `lrmodel` class
    * Use `modelCLR` to capture the clr-transformation rule
    * Use `predict` to deploy this rule to new data
* Modified backend code
    * Added `corRcpp` function from `correlateR` package
    * Added `linRcpp` function for Z-transformation
    * Added `lltRcpp` and `urtRcpp` to retrieve a half-matrix
    * Added `labRcpp` to label a half-matrix
* New `prop2prob` function
    * Allows hypothesis testing of rho equals naught
    * Tests differential proportionality

## propr 2.0.6
---------------------
* Modified visualization tools
    * `plotCheck` extended to all plot functions
    * `plot` method now calls `smear` function
    * `dendrogram` plot now rendered using `ggplot2`
    * `snapshot` plot now rendered using `ggplot2`
    * `bokeh` plot now on positive log scale
    * `plotly` support added
* Modified backend code
    * Temporarily removed `a_bool` function
* Modified `[` method
    * Removed `bool` and `copy` arguments

## propr 2.0.5
---------------------
* Modified backend code
    * New `a_bool` function returns thresholded boolean matrix
* Modified `[` method
    * New `bool` argument toggles whether to use `a_bool`
    * New `tiny` argument toggles whether to use `simplify`
    * New `copy` argument toggles `a_bool` copy-on-modify

## propr 2.0.4
---------------------
* New visualization tools
    * `slate` returns a table of VLR, VLS, and rho
    * `bokeh` plots pairs by the individual variances
* Modified index-naive plot functions
    * Now uses `fastcluster::hclust` implementation
    * New `prompt` argument turns off big data prompt
    * `prism` now depends on `slate` function
* Modified `dendrogram` function
    * Now uses `fastcluster::hclust` implementation
    * Now returns an `hclust` object

## propr 2.0.3
---------------------
* Modified `subset` method
    * Argument `select` now correctly rearranges features
* Modified `rhoRcpp` function
    * Now accommodates new `perb` function feature
* Modified `perb` function
    * New `select` argument returns subsetted matrix
    * This subset does not alter values of rho

## propr 2.0.2
---------------------
* Modified `perb` function
    * User can now specify name of `ivar` reference
* Altered `image` method
    * Now includes dendrogram with heatmap
    * No longer uses index pairs
    * Now called `snapshot`
* New `prism` function
* New `bucket` function
* New `mds` function
* New vignette

## propr 2.0.1
---------------------
* Modified `phit`, `perb` functions
    * These functions now force zero removal
* New `simplify` function
    * Subsets `propr` object based on index in `@pairs` slot
    * Returns an updated index

## propr 2.0.0
---------------------
* Modified `phit`, `perb` functions
    * Permutation testing removed
    * Added `lazyPairs` construct
        * Slot @pairs not populated until after `[`
* Modified `propr` Class
    * `@pairs` slot now integer vector
        * Populated with `indexPairs` function
        * Translate with `indexToCoord` function
    * `show` method updated for `lazyPairs` construct
    * `[` method completely redesigned
        * First argument specifies operation
        * Second argument specifies reference
        * Indexes @matrix based on these
    * `subset` method revised but still copy-on-modify
        * Resets @pairs when called
    * `$` method removed
* Visualization tools revised
    * `plot`, `image`, `dendrogram` methods
        * Improved performance
        * Compatible with new @pairs indexing
        * No longer requires column names
* Modified backend code
    * Rephrased code for `proprPhit`
    * Rephrased code for `proprPerb`
    * Rephrased code for `proprVLR`
    * All functions translated into C++
        * Estimated 80% reduction in RAM overhead
        * Estimated 100-fold performance increase
        * ALR methods no longer drop dimension
        * All have modify-in-place behavior

## propr 1.1.0
---------------------
* New orientation expected for input data
    * Updated backend and vignette accordingly
    * Removed redundant transpositions
* Fixed rare subsetting errors
* Tweaked plot methods

## propr 1.0.0
---------------------
* Introduced `phit` function
    * Implements Lovell's phi proportionality metric
    * Returns object of class `propr`
* Introduced `perb` function
    * Implements Erb's rho proportionality metric
    * Returns object of class `propr`
* Introduced `propr` Class
    * `show` method
        * Subsets `propr` based on `@pairs` slot
    * `subset` method
        * Subsets `propr` based on `@matrix` slot
    * `plot` method
        * Plots pairwise lr proportionality
    * `dendrogram` method
        * Plots clusters of lr-transformed data
    * `image` method
        * Plots heatmap of lr-transformed data
