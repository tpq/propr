#' @param counts A data matrix representing counts.
#'   It is assumed that the matrix contains numerical values only.
#' @param group A character vector representing group labels indicating the
#'  assignment of each count to different groups.
#' @param alpha The alpha parameter used in the alpha log-ratio transformation.
#' @param p The number of permutations to perform for calculating the false
#'  discovery rate (FDR). The default is 0.
#' @param weighted A logical value indicating whether weighted calculations
#'  should be performed. 
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
                  shrink = FALSE) {
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

  # Throw error if scenario not supported
  if (shrink && weighted) {
    stop("Shrinkage is not available for weighted computation yet.")
  }

  # Special handling for equivalent args
  if (identical(alpha, 0))
    alpha <- NA

  # Initialize @active, @weighted
  result <- new("propd")
  result@active <- "theta_d" # set theta_d active by default
  result@weighted <- weighted
  result@shrink <- shrink
  result@dfz <- 0

  # Initialize @counts, @group, @alpha
  result@counts <- as.data.frame(counts)
  result@group <- as.character(group)
  result@alpha <- as.numeric(alpha)
  result@permutes <- data.frame()

  ##############################################################################
  ### CALCULATE THETA WITH OR WITHOUT WEIGHTS
  ##############################################################################

  # Initialize @results
  result@results <-
    calculate_theta(
      result@counts,
      result@group,
      result@alpha,
      weighted = weighted,
      shrink = shrink
    )
  result@results$Zeros <- ctzRcpp(counts) # count number of zeros
  result@results$theta <-
    round(result@results$theta, 14) # round floats to 1

  # permute data
  if (p > 0) result <- updatePermutes(result, p)

  ##############################################################################
  ### GIVE HELPFUL MESSAGES TO USER
  ##############################################################################

  message("Alert: Use 'setActive' to select a theta type.")
  message("Alert: Use 'updateCutoffs' to calculate FDR.")
  message("Alert: Use 'updateF' to calculate F-stat.")

  return(result)
}