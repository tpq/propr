#' @param counts A data matrix representing counts.
#'   It is assumed that the matrix contains numerical values only.
#' @param group A character vector representing group labels indicating the
#'  assignment of each count to different groups.
#' @param alpha The alpha parameter used in the alpha log-ratio transformation.
#' @param p The number of permutations to perform for calculating the false
#'  discovery rate (FDR). The default is 0.
#' @param weighted A logical value indicating whether weighted calculations
#'  should be performed. If \code{TRUE}, the function will use limma-based
#'  weights for the calculations.
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
                  weighted = FALSE) {
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

  # Special handling for equivalent args
  if (identical(alpha, 0))
    alpha <- NA

  ##############################################################################
  ### OPTIONALLY REPLACE ZEROS AND SET UP propd OBJECT
  ##############################################################################

  if (is.na(alpha)) {
    ct <- simple_zero_replacement(counts)
  } else{
    ct <- counts
  }

  # Initialize @active, @weighted
  result <- new("propd")
  result@active <- "theta_d" # set theta_d active by default
  result@weights <- as.matrix(NA)
  result@weighted <- weighted
  result@dfz <- 0

  # Initialize @counts, @group, @alpha
  result@counts <- as.data.frame(ct)
  result@group <- as.character(group)
  result@alpha <- as.numeric(alpha)
  result@permutes <- data.frame()

  ##############################################################################
  ### CALCULATE THETA WITH OR WITHOUT WEIGHTS
  ##############################################################################

  # Initialize @weights
  if (weighted) {
    message("Alert: Calculating limma-based weights.")
    packageCheck("limma")
    design <-
      stats::model.matrix(~ . + 0, data = as.data.frame(group))
    v <- limma::voom(t(counts), design = design)
    result@weights <- t(v$weights)
  }

  # Initialize @results
  result@results <-
    calculate_theta(
      result@counts,
      result@group,
      result@alpha,
      weighted = result@weighted,
      weights = result@weights
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
