#' Get Matrix from Object
#'
#' This function provides a unified wrapper to retrieve a matrix
#'  of \code{propr} or \code{propd} values.
#'
#' @inheritParams getResults
#'
#' @return A matrix.
#'
#' @export
getMatrix <- function(object) {
  NVTX_PUSH("getMatrix", 0)

  if (class(object) == "propr") {
    NVTX_PUSH("getMatrix_propr", 0)
    mat <- object@matrix
    NVTX_POP()

  } else if (class(object) == "propd") {
    NVTX_PUSH("getMatrix_propd", 0)
    mat <- results_to_matrix(object@results, features = colnames(object@counts))
    NVTX_POP()

  } else {
    NVTX_POP()  # pop "getMatrix" before error
    stop("Provided 'object' not recognized.")
  }

  NVTX_POP()    # pop "getMatrix"
  return(mat)
}

#' Get Matrix from Results
#' 
#' This function converts the results data frame into a matrix.
#' 
#' @param results A \code{data.frame} of results.
#' @param what A character string. The column name of the results data frame to be converted into a matrix.
#' @param features A vector of features. Default is NULL.
#' 
#' @return A matrix.
#' 
#' @export
results_to_matrix <- function(results, what = 'theta', features = NULL) {
  NVTX_PUSH("results_to_matrix", 0)

  NVTX_PUSH("determine_features_and_indices", 0)
  # if pair and partner are already named
  if (!is.numeric(results$Pair) && !is.numeric(results$Partner)) {
    if (is.null(features)) {
      features <- unique(c(results$Pair, results$Partner))
    }
    nfeatures <- length(features)
    pair <- match(results$Pair, features)
    partner <- match(results$Partner, features)
    if (any(is.na(pair)) || any(is.na(partner))) {
      NVTX_POP()  # determine_features_and_indices
      NVTX_POP()  # results_to_matrix
      stop("Some features are not found in the results data frame.")
    }

  # if pair and partner are still indices
  } else {
    if (is.null(features)) {
      features <- sort(unique(c(results$Pair, results$Partner)))
      nfeatures <- max(features)
    } else {
      if (length(features) != max(results$Pair, results$Partner)) {
        NVTX_POP()  # determine_features_and_indices
        NVTX_POP()  # results_to_matrix
        stop("The length of 'features' does not match the number of features in the results data frame.")
      }
      nfeatures <- length(features)
    }
    pair <- results$Pair
    partner <- results$Partner
  }
  NVTX_POP()  # determine_features_and_indices

  NVTX_PUSH("build_matrix", 0)
  # convert the results data frame into a matrix
  mat <- vector2mat(results[, what], pair, partner, nfeatures)
  diag(mat) <- 0
  if (!is.numeric(features)) rownames(mat) <- colnames(mat) <- features
  NVTX_POP()  # build_matrix

  NVTX_POP()  # results_to_matrix
  return(mat)
}
