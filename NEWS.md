## propr 3.1.8
---------------------
* Modified `propd` methods
    * Users can now disable alpha transformation by setting `alpha = 0`
    * Improve `NaN` theta value replacement with 1
    * Modified `updateF` for `moderated = TRUE`
        * Now offsets counts by 1 to prevent zeros in reference set
        * Now correctly checks for zeros in reference set

## propr 3.1.7
---------------------
* Check `propr` and `propd` input for negative counts
* Modified `propd` methods
    * Add `dfz` slot to the `propd` class (defaults to 0)
        * Have `updateF` populate `dfz` slot if moderated
        * Use `dfz` to calculate p-value from F-stat
    * Implement revised F-stat moderation
        * Now uses simple moderation term based only on LRV
        * Update `test-Fstat.R` to reflect changes

## propr 3.1.6
---------------------
* Modified `propd` methods
    * Allow `p = 0` when initializing `propd` objects
    * Have `calculateTheta` compute new weights each permutation

## propr 3.1.5
---------------------
* Update `propr` to work with Rcpp >= 0.12.12
* Modified `propd` methods
    * All out-of-bounds theta_mod replaced with 1
    * Add `@Fivar` slot to the `propd` class
    * Extend `updateCutoffs` to `theta_mod`

## propr 3.1.4
---------------------
* Automatically set `@matrix` column names when using `propr`
* Address notice by CRAN about using Suggests conditionally
    * Vignette will no longer calculate theta when `weighted = TRUE`
    * Tests for "Fstat" and "theta" now check for `limma`

## propr 3.1.3
---------------------
* New `ratios` and `ratiosRcpp` functions recast matrix as feature ratios

## propr 3.1.2
---------------------
* Rename "Calculating Differential Proportionality" vignette
* Update CITATION file and README
* Modified `propd` methods
    * Stop plot functions from creating surplus figure in notebook
    * Use theta_d network method for any non-emergent theta type
    * Add `suppressWarnings` to `compositions::plot` call

## propr 3.1.1
---------------------
* New vignette to discuss F-statistics
* Modified `propd` methods
    * Add `suppressWarnings` to `compositions::acomp` call

## propr 3.1.0
---------------------
* Modified `propd` methods
    * The `propd` function no longer calls `updateCutoffs`
    * The `calculateTheta` function no longer calculates F-statistic
    * Update vignette to reflect these changes
    * Object now stores weights if weighted
* Modified `calculateTheta`
    * Added `weights` argument to pass pre-computed weights
    * Still calculates weights if `NA` weights argument
    * Consolidated weighted and alpha lrv calls
* Modified C++ backend
    * Extend `lrv` to weighted alpha calls
    * Added unit test for weighted alpha calls
    * Merge `boxRcpp` with `lrv` and remove `boxRcpp`
    * Merge `lrmRcpp` with `lrm` and remove `lrmRcpp`
    * Rename `lrvMod` function to `omega`
* New `updateF` function
    * Added new `ivar2index` function used by `propr` and `updateF`
    * Calculates 4 types of moderated F-statistics
    * Calculates 4 types of non-moderated F-statistics
    * Update `propd` documentation to reflect change
    * Extend moderation to non-clr references

## propr 3.0.7
---------------------
* New `corr` function calculates log-ratio based correlation
* Fixed mistake with vignette indexing

## propr 3.0.6
---------------------
* Modified visualization tools
    * Correctly spell `geiser` function as `geyser`
    * Allow `propd` methods to accept 0s
* Update documentation
    * Finish revision of all vignettes

## propr 3.0.5
---------------------
* Modified `propd` Class
    * Now only replaces 0s if alpha is missing or `NA`
    * Now tabulates 0s whether replaced or not
* Modified `calculateTheta`
    * Now replaces 0s to calculate LRM if alpha is provided
    * Now alerts user when using alpha to approximate LRV
* Update documentation
    * Create new DESCRIPTION that includes differential proportionality
    * Create new README that includes differential proportionality
    * Revise "a_introduction" and move portion to stand-alone critique
    * Revise "b_visualization" and fix co-cluster selection mistake
    * Fork off and revise "d_advanced" as stand-alone critique
    * Spell check and revise other vignettes

## propr 3.0.4
---------------------
* Modified `calculateTheta`
    * Now saves log-ratio variance (LRV) modifier in output
* Modified visualization tools
    * New `decomposed` function for LRV decomposition
    * Uses LRV modifier for weighted theta types

## propr 3.0.3
---------------------
* Modified visualization tools
    * Updated `plot.propd` method to display `theta_f`

## propr 3.0.2
---------------------
* Modified `propd` Class
    * Added `@weighted` slot now used by `updateCutoffs`
* Implement # Initialize @theta for "weighted theta" calculation
    * Added `wtmRcpp` for weighted mean calculation
    * Added `wtvRcpp` for weighted variance calculation
    * Added `lrm` function with optional weighted calculations
    * Added `lrv` function with optional weighted calculations
* Implement "weighted theta" calculation
    * Added `calculateThetaW_old` for unit testing
    * Added weighted lrv calculation
* Modified `calculateTheta`
    * Moved log-ratio mean calculations here
    * Added weighted lrm calculation

## propr 3.0.1
---------------------
* Modified `propd` Class
    * Added `@active` slot now used by `updateCutoffs`
    * Added `setActive` method to switch between theta types
* Modified `calculateTheta`
    * Added `theta_f` which equals `1 - theta_e`
    * Added `only` argument to retrieve only one theta type

## propr 3.0.0
---------------------
* Modified package skeleton
    * Added differential proportionality article to CITATION
    * Added propriety project pre-merger changelog to OLDNEWS.md
* Manually merged R functions from propriety project
    * New functions: `propd` and methods, `calculateTheta`, `updateCutoffs`
    * Added deprecated functions anticipating unit tests
* Manually merged C++ functions from propriety project
    * New functions: `half2matrix`, `boxRcpp`, `ctzRcpp`, `lrmRcpp`
    * Added deprecated functions anticipating unit tests
* Modified unit tests
    * Added `requireNamespace` check for `ALDEx2` tests
    * Added unit tests from propriety project
* Modified data
    * Rebuilt `pd.d` and `pd.e` data
    * Added `top.counts` as count data using filtered `caneToad.counts`
    * Added `top.lr` as log-ratio data using filtered `caneToad.counts`
    * Removed `top` data

## propr 2.2.0
---------------------
* Modified dataset from Rollins et al. 2015
    * Removed transcripts with < 10 counts in < 10 samples
* Add dataset from Marguerat et al. 2012
    * Absolute data stored as `marg.abs`

## propr 2.1.9
---------------------
* Replace `rhoToPhs` function with `lr2phs`
    * Fixed `phis` bug when using alr-transformation
    * Fixed `aldex2propr` bug when using alr-transformation
    * Add unit tests for `lr2` with `ivar`
* Create `initialize` method to handle `ivar` and `select`
    * `phit`, `perb`, and `phis` use `new` with `lr2`
    * `phit` now accepts `ivar` and `select`
    * `phit` now has `pca` and `snapshot`
    * Add `iqlr` and "multi"-alr
* Updated "Frequently Asked Questions" vignette

## propr 2.1.8
---------------------
* Modified visualization tools
    * The `bucket` function now uses non-parametric ANOVA
    * The `cytescape` function now supports 3D visualization
    * The `cytescape` function now names columns correctly
* Add "Frequently Asked Questions" vignette
* Removed `lrmodel` module

## propr 2.1.7
---------------------
* Force column and row names when calculating proportionality
    * Added `plotCheck` to `cytescape` function

## propr 2.1.6
---------------------
* The `aldex.cor` function now returns average p-value
    * Now uses `cor.test` as `lr2cor` foundation instead of `cor`
* Modified `aldex2propr` function
    * Added `select` argument like in `perb` function
    * Documented `select` argument in Details
* Modified `perb` function
    * Alerts user when 'ivar' is missing from 'select'
    * Documented `select` argument in Details

## propr 2.1.5
---------------------
* Stop exporting `progress` and `migraph` functions
* Color `cytescape` edges correctly when plotting phi or phs
* Add `progress` bar to `smear` function

## propr 2.1.4
---------------------
* Add `progress` bar to `aldex2propr` and `aldex.cor` functions
* Penalize `lr2cor` and `aldex.cor` p-values for two-tailed test
* Add the `migraph` module to help make `igraph` networks
* Modified visualization tools
    * Rebuilt `cytescape` using `migraph` module

## propr 2.1.3
---------------------
* New `aldex2propr` function converts`aldex.clr` object
    * Uses Monte Carlo instances from `ALDEx2` package
* New `lr2cor` and `aldex.cor` functions
    * Measures feature associations with continuous variables
* Modified backend code
    * New `backend.h` allows import of `backend.cpp`
    * New `lr2` functions calculate proportionality
        * Input log-ratio transformed counts
* Modified visualization tools
    * Removed `minPairs` argument from `cytescape`

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
