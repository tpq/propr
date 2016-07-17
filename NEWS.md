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
  * Implements Lovell's \phi proportionality metric
  * Returns object of class `propr`
* Introduced `perb` function
  * Implements Erb's \rho proportionality metric
  * Returns object of class `propr`
* Introduced `propr` Class
  * `show` method
    * Subsets `propr` based on `@pairs` slot
  * `subset` method
    * Subsets `propr` based on `@matrix` slot
  * `plot` method
    * Plots pairwise *lr proportionality
  * `dendrogram` method
    * Plots clusters of *lr-transformed data
  * `image` method
    * Plots heatmap of *lr-transformed data
