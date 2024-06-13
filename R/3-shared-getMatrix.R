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
    df <- object@matrix
  }else if(class(object) == "propd"){
    df <- get_theta_matrix(object)
  }else{
    stop("Provided 'object' not recognized.")
  }

  return(df)
  }
