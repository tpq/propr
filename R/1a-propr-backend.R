#' Index Reference for Selecting Features
#'
#' This function computes an index reference for selecting features from the
#'  input count matrix based on the provided `ivar` argument. The index
#'  reference can be obtained using different options, such as selecting all
#'  features, using a user-defined vector of indices, selecting features by
#'  name, using the "clr" transformation, or computing the IQLR
#'  (Interquartile Range Log Ratio) features.
#'
#' @inheritParams propr
#' @return A numeric vector representing the indices of features selected.
#' @examples
#' # Sample input count data
#' data <- matrix(c(1, 2, 3, 4, 0, 6), nrow = 2, byrow = TRUE)
#'
#' # Index reference using all features
#' all_features <- index_reference(data, ivar = "all")
#'
#' # Index reference using custom numeric vector
#' custom_indices <- index_reference(data, ivar = c(1, 3))
#'
#' # Index reference using IQLR features
#' iqlr_features <- index_reference(data, ivar = "iqlr")
#'
#' @export
index_reference <- function(counts, ivar) {
  `%is%` <- function(a, b)
    identical(a, b)

  if (any(is.na(ivar))) {
    stop("DEBUG IVAR001: NA 'ivar' should not appear here.")
  } else if (ivar %is% 0 |
             ivar %is% NULL | ivar %is% "all" | ivar %is% "clr") {
    use <- 1:ncol(counts) # use all features for geometric mean

  } else if (ivar %is% "iqlr") {
    use <- compute_iqlr(counts)

  } else if (is.character(ivar)) {
    use <-
      which(colnames(counts) %in% ivar) # use features by name

  } else{
    use <- sort(ivar) # use features given by number
  }

  return(use)
}

#' Compute IQLR (Interquartile Range Log Ratio) Features
#'
#' This function computes the IQLR features from the input count matrix.
#'  The IQLR is based on the log-ratio transformation of the counts, and
#'  it selects features with variance values falling within the
#'  interquartile range.
#'
#' @inheritParams propr
#' @return A numeric vector representing the indices of features selected by IQLR.
#' @examples
#' # Sample input count data
#' data <- matrix(c(1, 2, 3, 4, 0, 6), nrow = 2, byrow = TRUE)
#'
#' # Compute IQLR features
#' iqlr_features <- compute_iqlr(data)
#'
#' @export
compute_iqlr <- function(counts) {
  counts <- suppressMessages(simple_zero_replacement(counts))
  counts.clr <- apply(log(counts), 1, function(x) {
    x - mean(x)
  })
  counts.var <- apply(counts.clr, 1, stats::var)
  quart <-
    stats::quantile(counts.var) # use features with unextreme variance
  use <- which(counts.var < quart[4] & counts.var > quart[2])
  return(use)
}

#' Simple Zero Replacement in a Count Matrix
#'
#' This function replaces zeros with the next smallest non-zero value in the
#'  input count matrix. If the matrix contains no zeros, it produces an
#'  informational message indicating that no replacements were made.
#'
#' @inheritParams logratio_with_alpha
#' @return A matrix with zero values replaced by the next smallest non-zero value.
#' @examples
#' # Sample input count data with zeros
#' data <- matrix(c(0, 2, 3, 4, 5, 0), nrow = 2, byrow = TRUE)
#'
#' @export
simple_zero_replacement <- function(ct) {
  if (any(ct == 0)) {
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])
    message("Alert: replacing zeros with minimun value.")
  } else{
    message("Alert: No 0s found that need replacement.")
  }

  return(ct)
}

#' Perform log-ratio transformation without alpha parameter
#'
#' This function applies a log-ratio transformation to a given data matrix
#'  without using an alpha parameter. The log-ratio transformation is based
#'  on a selected subset of columns specified by the `use` argument.
#'
#' @inheritParams logratio_with_alpha
#' @return A matrix containing the log-ratio transformed data.
#' @examples
#' # Sample input data
#' data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#'
#' # Applying log-ratio transformation to rows using columns 2 and 3
#' result <- logratio_without_alpha(data, use = c(2, 3))
#'
#' @export
logratio_without_alpha <- function(ct, use) {
  # Use g(x) = Mean[log(x)] to log-ratio transform data
  logX <- log(ct)
  logSet <- logX[, use, drop = FALSE]
  ref <- rowMeans(logSet)
  lr <- sweep(logX, 1, ref, "-")
  return(lr)
}

#' Perform log-ratio transformation with alpha parameter
#'
#' This function applies a log-ratio transformation to a given data matrix
#'  using an alpha parameter. The log-ratio transformation is based on a
#'  selected subset of columns specified by the `use` argument.
#'  The transformation formula is: \code{log(x/ref) = log(x) - log(ref) =
#'  [x^alpha-1]/alpha - [ref^alpha-1]/alpha}, where \code{x} represents the
#'  data matrix and \code{ref} is the reference value calculated as the mean
#'  of the selected subset of columns.
#'
#' @param ct A data matrix for which the log-ratio transformation will be
#'  performed. It is assumed that the matrix contains numerical values only.
#' @param use An integer vector specifying the subset of columns to be used
#'  for the log-ratio transformation.
#' @param alpha The alpha parameter used in the transformation.
#'
#' @return A matrix containing the log-ratio transformed data.
#'
#' @examples
#' # Sample input data
#' data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#'
#' # Applying log-ratio transformation with alpha = 2 to rows using columns 2 and 3
#' result <- logratio_with_alpha(data, use = c(2, 3), alpha = 2)
#'
#' @export
logratio_with_alpha <- function(ct, use, alpha) {
  # Since y = log(x) = [x^a-1]/a, ref = Mean[log(x)] = Mean[y]
  # Calculate log(x/ref) = log(x) - log(ref) = [x^a-1]/a - [ref^a-1]/a
  aX <- (ct ^ alpha - 1) / alpha
  aSet <- aX[, use, drop = FALSE]
  ref <- rowMeans(aSet)
  lr <- sweep(aX, 1, ref, "-")
  return(lr)
}

#' Perform log-ratio transformation
#'
#' This function applies a log-ratio transformation to a given data matrix
#'  with or without using an alpha parameter.
#'
#' @inheritParams propr
#' @return A matrix with log-ratio transformed values.
#' @examples
#' # Sample counts matrix
#' counts_matrix <- matrix(c(10, 20, 30, 40, 0, 50, 60, 70, 0), nrow = 3, byrow = TRUE)
#' colnames(counts_matrix) <- c("A", "B", "C")
#' rownames(counts_matrix) <- c("Sample1", "Sample2", "Sample3")
#'
#' # Perform log-ratio transformation without alpha
#' logratio_without_alpha(counts_matrix, use = "A")
#'
#' # Perform log-ratio transformation with alpha 0.5
#' logratio_with_alpha(counts_matrix, use = "A", alpha = 0.5)
#'
#' # Skip log-ratio transformation
#' logratio(counts_matrix, ivar = NA)
#'
#' @export
logratio <- function(counts, ivar, alpha) {
  counts <- as_safe_matrix(counts)
  if (is.na(ivar)) {
    message("Alert: Skipping built-in log-ratio transformation.")
    ct <- simple_zero_replacement(counts)
    lr <- ct
  } else {
    use <- index_reference(counts, ivar) # Index reference
    if (is.na(alpha)) {
      ct <- simple_zero_replacement(counts)
      lr <- logratio_without_alpha(ct, use)
    } else {
      lr <- logratio_with_alpha(counts, use, alpha)
    }
  }
  return(lr)
}

#' Covariance Shrinkage and Partial Correlation Calculation
#'
#' This function performs covariance shrinkage and calculates the partial
#'  correlation matrix based on the input count data matrix. The function can
#'  output the results in two formats: centered log-ratio (clr) or
#'  additive log-ratio (alr).
#'
#' @param M A data matrix representing the counts. Each row should represent
#'  observations, and each column should represent different variables.
#' @param outtype A character vector specifying the output type. It can take
#'  either "clr" (centered log-ratio) or "alr" (additive log-ratio).
#'
#' @return A matrix representing the partial correlation matrix.
#'
#' @examples
#' # Sample input count data
#' data <- iris[,1:4]
#'
#' # Calculate partial correlation matrix using clr transformation
#' result_clr <- basis_shrinkage(data, outtype = "clr")
#'
#' # Calculate partial correlation matrix using alr transformation
#' result_alr <- basis_shrinkage(data, outtype = "alr")
#'
#' @export
basis_shrinkage <- function(M, outtype = c("clr", "alr")) {
  packageCheck("corpcor")
  outtype <- match.arg(outtype)

  # transform counts to log proportions
  P <- M / rowSums(M)
  B <- log(P)

  # covariance shrinkage
  D  <- ncol(M)
  Cb <- corpcor::cov.shrink(B, verbose = FALSE)
  if (outtype == "alr") {
    F   <- cbind(diag(rep(1, D - 1)), rep(-1, D - 1))
    Cov <- F %*% Cb %*% t(F)
  } else if (outtype == "clr") {
    G   <- diag(rep(1, D)) - matrix(1 / D, D, D)
    Cov <- G %*% Cb %*% G
  }

  # partial correlation
  PC <- corpcor::cor2pcor(Cov)

  # make output to have same dimensions as input
  if (outtype == "alr") {
    PC <- cbind(PC, replicate(nrow(PC), 0))
    PC <- rbind(PC, replicate(ncol(PC), 0))
    PC[nrow(PC), ncol(PC)] <- 1
  }

  return(PC)
}

#' Wrapper to check if a metric is 'direct' or 'positive'
#' 
#' If a metric is better when increased, it is considered 'direct'.
#' This function checks if a propr metric is direct.
#' 
#' @param metric A character vector specifying the metric to check.
#' @return A logical value indicating whether the metric is direct.
#' @noRd
metric_is_direct <- function(metric) {
    metrics <- c("rho", "cor", "pcor", "pcor.shrink", "pcor.bshrink")
    return(metric %in% metrics)
  }