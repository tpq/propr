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

  if(class(object) == "propr"){
    mat <- object@matrix

  }else if(class(object) == "propd"){
    mat <- results_to_matrix(object@results, features=colnames(object@counts))

  }else{
    stop("Provided 'object' not recognized.")
  }

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
results_to_matrix <- function(results, what='theta', features = NULL) {

  # if pair and partner are already named
  if (!is.numeric(results$Pair) && !is.numeric(results$Partner)) {
    if (is.null(features)) {
      features <- unique(c(results$Pair, results$Partner))
    }
    nfeatures <- length(features)
    pair <- match(results$Pair, features)
    partner <- match(results$Partner, features)
    if (any(is.na(pair)) || any(is.na(partner))) {
      stop("Some features are not found in the results data frame.")
    }

  } else {
    if (is.null(features)) {
      features <- sort(unique(c(results$Pair, results$Partner)))
      nfeatures <- max(features)
    } else {
      if (length(features) != max(results$Pair, results$Partner)) {
        stop("The length of 'features' does not match the number of features in the results data frame.")
      }
      nfeatures <- length(features)
    }
    pair <- results$Pair
    partner <- results$Partner
  }

  # convert the results data frame into a matrix
  mat <- vector2mat(results[,what], pair, partner, nfeatures)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- features

  return(mat)
}