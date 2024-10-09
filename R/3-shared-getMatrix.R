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
    mat <- results_to_matrix(object@results)

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
#' 
#' @return A matrix.
#' 
#' @export
results_to_matrix <- function(results, what='theta') {

  # if the results data frame has the Pair and Partner columns as names
  if (!is.numeric(results$Partner)) {
    features <- unique(c(results$Pair, results$Partner))
    pair <- match(results$Pair, features)
    partner <- match(results$Partner, features)
    nfeatures <- length(features)

  # if the results data frame has the Pair and Partner columns as indices
  } else {
    features <- sort(unique(c(results$Pair, results$Partner)))
    features <- colnames(pr@counts)[features]
    pair <- results$Pair
    partner <- results$Partner
    nfeatures <- length(features)
  }

  # convert the results data frame into a matrix
  mat <- vector2mat(results[,what], pair, partner, nfeatures)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- features

  return(mat)
}