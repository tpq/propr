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

  if (length(use) == 0) {
      stop("No reference detected. Please provide a valid 'ivar' argument.")
  } else if(length(use) > ncol(counts)) {
    stop("Detected more reference variables than existing in the data. Please provide a valid 'ivar' argument.")
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
#' # Applying log-ratio transformation to rows using columns 2 and 3 as reference
#' result <- logratio_without_alpha(data, use = c(2, 3))
#'
#' @export
logratio_without_alpha <- function(ct, use) {
  if (any(ct == 0)){
    stop("Error: please first replace the zeros before running logratio.")
  }
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
#' counts_matrix <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80, 90), nrow = 3, byrow = TRUE)
#' colnames(counts_matrix) <- c("A", "B", "C")
#' rownames(counts_matrix) <- c("Sample1", "Sample2", "Sample3")
#'
#' # Perform log-ratio transformation without alpha
#' logratio(counts_matrix, ivar = "A")
#'
#' # Perform log-ratio transformation with alpha 0.5
#' logratio(counts_matrix, ivar = "A", alpha = 0.5)
#'
#' # Skip log-ratio transformation
#' logratio(counts_matrix, ivar = NA)
#'
#' @export
logratio <- function(counts, ivar, alpha=NA) {
  counts <- as_safe_matrix(counts)
  if (length(ivar) == 1 && is.na(ivar)) {
    message("Alert: Skipping built-in log-ratio transformation.")
    lr <- counts
  } else {
    use <- index_reference(counts, ivar) # Index reference
    if (is.na(alpha)) {
      lr <- logratio_without_alpha(counts, use)
    } else {
      lr <- logratio_with_alpha(counts, use, alpha)
    }
  }
  return(lr)
}

#' Basis Covariance Shrinkage and Partial Correlation Calculation
#'
#' This function performs covariance shrinkage on the basis matrix
#'  and calculates the partial correlation matrix. The function can
#'  output the results in two formats: centered log-ratio (clr) or
#'  additive log-ratio (alr).
#'
#' @param ct A data matrix representing the counts. Each row should represent
#'  observations, and each column should represent different variables.
#' @param outtype A character vector specifying the output type. It can take
#'  either "clr" (centered log-ratio) or "alr" (additive log-ratio). Since
#'  the reference does not affect the logratio partial correlation coefficients,
#'  the index of the reference is not needed for "alr". Also, "clr" is recommended
#'  because of the same reason, at the same time avoiding losing one dimension.
#'
#' @return A matrix representing the shrunk partial correlation matrix.
#'
#' @examples
#' # Sample input count data
#' data <- iris[,1:4]
#'
#' # Calculate partial correlation matrix using clr transformation
#' result_clr <- pcor.bshrink(data, outtype = "clr")
#'
#' # Calculate partial correlation matrix using alr transformation
#' result_alr <- pcor.bshrink(data, outtype = "alr")
#'
#' @export
pcor.bshrink <- function(ct, outtype = c("clr", "alr")) {
  NVTX_PUSH("pcor.bshrink", 0)
  packageCheck("corpcor")
  outtype <- match.arg(outtype)

  # transform counts to log proportions
  logP <- log( ct / rowSums(ct) )

  # covariance shrinkage
  covB <- corpcor::cov.shrink(logP, verbose = FALSE)
  lambda <- attr(covB, "lambda")

  # convert basis covariance matrix to clr/alr covariance matrix
  NVTX_PUSH("basis covariance matrix to clr/alr", 0)
  D  <- ncol(ct)
  if (outtype == "alr") {
    F   <- cbind(diag(rep(1, D - 1)), rep(-1, D - 1))
    cov <- F %*% covB %*% t(F)
  } else if (outtype == "clr") {
    G   <- diag(rep(1, D)) - matrix(1 / D, D, D)
    cov <- G %*% covB %*% G
  }
  NVTX_POP()

  # partial correlation
  NVTX_PUSH("corpcor::cor2pcor", 0)
  pcor <- corpcor::cor2pcor(cov)
  NVTX_POP()

  # make output to have same dimensions as input
  # alr partial correlation has one less dimension,
  # so here we add a row and a column of 0s
  if (outtype == "alr") {
    pcor <- cbind(pcor, 0)
    pcor <- rbind(pcor, 0)
    pcor[ncol(pcor), ncol(pcor)] <- 1
  }
  NVTX_POP()
  return(list(matrix = pcor, lambda = lambda))
}