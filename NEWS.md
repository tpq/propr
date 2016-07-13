## propr 1.1.0.9000
---------------------
* Added `lazyPairs` argument to `phit`
  * Does not populate @pairs slot until after subset
* Added `lazyPairs` argument to `perb`
  * Does not populate @pairs slot until after subset
* Adjusted `show` for `lazyPairs` update
//* Adjusted `[` and `subset` for `lazyPairs` update
//* Adjusted `$` for `lazyPairs` update
//* Adjusted `plot` for `lazyPairs` update
* Rephrased code for `proprPhit`
* Rephrased code for `proprPerb`
* Rephrased code for `proprVLR`

## propr 1.1.0
---------------------
* New orientation expected for input data
  * Updated back-end and vignette accordingly
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
