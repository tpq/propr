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
getMatrix <-
  function(object) {
    df <- object@matrix
    return(df)
  }
