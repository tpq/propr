#' @param counts A data matrix representing counts.
#'   It is assumed that the matrix contains numerical values only.
#' @param group A character vector representing group labels indicating the
#'  assignment of each count to different groups.
#' @param alpha The alpha parameter used in the alpha log-ratio transformation.
#' @param p The number of permutations to perform for calculating the false
#'  discovery rate (FDR). The default is 0.
#' @param weighted A logical value indicating whether weighted calculations
#'  should be performed. 
#' @param weights A custom matrix of weights.
#' @param shrink A logical value indicating whether to apply shrinkage
#' 
#' @return A \code{propd} object containing the computed theta values,
#'  associated count matrix, group labels, and other calculated statistics.
#'
#' @details The \code{propd} function creates a \code{propd} object, which is
#'  used for differential analysis of regulatory relationships between features.
#'  It performs log-ratio transformation, calculates the variance of log-ratio
#'  (VLR), and computes the theta values using different formulas. The object
#'  stores the count matrix, group labels, alpha parameter, and other relevant
#'  information needed for further analysis and visualization.
#'
#' @examples
#' # Sample input count data and group assignments
#' data <- iris[1:100, 1:4]
#' group <- iris[1:100, 5]
#'
#' # Create a propd object for differential analysis
#' result <- propd(data, group, alpha = 0.5)
#'
#' @rdname propd
#' @export
propd <- function(counts,
                  group,
                  alpha = NA,
                  p = 0,
                  weighted = FALSE,
                  weights = as.matrix(NA),
                  shrink = FALSE) {
  nvtxR::nvtx_push_range("propd", 0)
  ##############################################################################
  ### CLEAN UP ARGS
  ##############################################################################

  # Clean "count matrix"
  counts <- as_safe_matrix(counts)

  # Clean group
  if (inherits(group, "factor"))
    group <- as.character(group)
  if (!inherits(group, "character"))
    stop("Provide group labels as a character vector.")
  if (length(group) != nrow(counts))
    stop("Too many or too few group labels.")

  # set weighted to TRUE if weights are provided
  if (!is.na(weights[1,1])) {
    weighted <- TRUE
    # error if permutation is requested
    if (p > 0) {
      stop("Permutation is not available with custom weights yet.")
    }
  }

  if (shrink && weighted) {
    stop("Shrinkage is not available for weighted computation yet.")
  }

  # Special handling for equivalent args
  if (identical(alpha, 0))
    alpha <- NA

  ##############################################################################
  ### OPTIONALLY REPLACE ZEROS AND SET UP propd OBJECT
  ##############################################################################

  if (is.na(alpha)) {
    nvtxR::nvtx_push_range("simple_zero_replacement", 1)
    ct <- simple_zero_replacement(counts)
    nvtxR::nvtx_pop_range()
  } else{
    ct <- counts
  }

  # Initialize @active, @weighted
  result <- new("propd")
  result@active <- "theta_d" # set theta_d active by default
  result@weighted <- weighted
  result@shrink <- shrink
  result@dfz <- 0

  # Initialize @counts, @group, @alpha
  result@counts <- as.data.frame(ct)
  result@group <- as.character(group)
  result@alpha <- as.numeric(alpha)
  result@permutes <- data.frame()

  ##############################################################################
  ### CALCULATE THETA WITH OR WITHOUT WEIGHTS
  ##############################################################################

  # Initialize @results
  nvtxR::nvtx_push_range("calculate_theta", 1)
  result@results <-
    calculate_theta(
      result@counts,
      result@group,
      result@alpha,
      weighted = weighted,
      weights = weights,
      shrink = shrink
    )
    nvtxR::nvtx_pop_range()
    
  nvtxR::nvtx_push_range("ctzRcpp", 1)
  result@results$Zeros <- ctzRcpp(counts) # count number of zeros
  nvtxR::nvtx_pop_range()

  result@results$theta <-round(result@results$theta, 14) # round floats to 1
  

  # permute data
  if (p > 0) {
      nvtxR::nvtx_push_range("updatePermutes", 1)
      result <- updatePermutes(result, p)
      nvtxR::nvtx_pop_range()
  }

  ##############################################################################
  ### GIVE HELPFUL MESSAGES TO USER
  ##############################################################################

  message("Alert: Use 'setActive' to select a theta type.")
  message("Alert: Use 'updateCutoffs' to calculate FDR.")
  message("Alert: Use 'updateF' to calculate F-stat.")
  nvtxR::nvtx_pop_range()
  return(result)
}
