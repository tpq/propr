## propriety 1.0.1
---------------------
* Replace `alphaTheta` with improved `calculateTheta`
    * Now calculates emergent proportionality (theta_e)
    * Now checks for `@theta$lrv == 0` (NaN theta)
* Modified `propd` function
    * FDR now handled by `updateCutoffs` function
    * New `setDisjointed` function makes theta_d active
    * New `setEmergent` function makes theta_e active
    * Now calculates F from theta_d
* Modified visualization tools
    * Remove "Bridged" and "Missing" pairs from figures
    * `geiser` function omits "Bridged" and "Missing" pairs
    * `gemini` function omits "Bridged" and "Missing" pairs
    * New `bowtie` function plots log-ratio means by group
* Modified `plot` function
    * New `plotSkip` argument used by `pals` to skip plot
    * Colors edges by LRM if using theta_d
    * Colors edges by LRV if using theta_e

## propriety 1.0.0
---------------------
* Include a sample `propd` object for vignette examples
    * Built using non-zero features with at least 40 counts
* Add conceptualization and implementation of PAL
    * New `pals` function computes the popular adjacent ligand
* Add table and visualization tools
    * New `shale` function produces table used for graphing
    * New `geiser` function plots VLR1 vs. VLR2
    * New `gemini` function plots log-fold VLR vs. LRMs
    * New `slice` function plots log-ratio abundance
* Modified `plot` function
    * Labels disjointedly proportional edges in yellow
    * Now augments graph using an indexed `propr` object
    * Integer cutoff now selects top N pairs
    * Now supports 3D visualization
* Add `propd` vignette

## propriety 0.0.7
---------------------
* Modified `calculateTheta` and `alphaTheta` functions
    * No longer returns a sorted `data.table` object

## propriety 0.0.6
---------------------
* Modified `calculateTheta` and `alphaTheta` functions
    * Now returns a sorted `data.table` object
* Add `progress` bar and `migraph` functions
* New `plot` method for `propd` object

## propriety 0.0.5
---------------------
* New `half2mat` function builds matrix from half-matrix
* New `propd2propr` function converts `propd` to `propr`
    * Allows `propd` to inherit `propr` methods

## propriety 0.0.4
---------------------
* Modified `propd` function
    * Use `ctzRcpp` to track joint zero frequency
    * Use `lrmRcpp` to calculate log-ratio mean
    * Help file includes messy LaTeX formulae
* Update documentation
    * README and vignette outlined with theory
* Deprecate outdated C++ code

## propriety 0.0.3
---------------------
* New `propd` class and function
    * Provides front-end for theta calculation
    * `alpha` argument toggles theta method
    * Include a unit test for new FDR
* Extended proportionality statistic theta
    * New `alphaTheta` function
    * Back-end handled by the `boxRcpp` function
* Included slow implementations for testing
    * Deprecated `alphaTheta_old`
* All permutation functions now deprecated

## propriety 0.0.2
---------------------
* Extended proportionality statistic theta
    * New `calculateFDR` function
* Removed unused C++ code

## propriety 0.0.1
---------------------
* Introduced proportionality statistic theta
    * New `calculateTheta` function
    * New `permuteTheta` function
* Included slow implementations for testing
    * Deprecated `calculateTheta_old`
    * Deprecated `permuteTheta_old`
* Added unit tests
    * New and old `calculateTheta` match
    * New and old `permuteTheta` match
