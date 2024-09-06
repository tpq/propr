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
    mat <- vector2mat(object@results$theta, object@results$Partner, object@results$Pair)
    rownames(mat) <- colnames(object@counts)
    colnames(mat) <- colnames(object@counts)
    diag(mat) <- 0

  }else{
    stop("Provided 'object' not recognized.")
  }

  return(mat)
  }
