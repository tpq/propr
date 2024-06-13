#' Get Object as Adjacency Matrix
#'
#' This function uses \code{getMatrix} to return a \code{propr}
#' or \code{propd} object as an adjacency matrix. The diagonal
#' is set to 0.
#'
#' @inheritParams getResults
#'
#' @return An adjacency matrix.
#'
#' @export
getAdj <- function(object, cutoff = 1, above = FALSE){
  m <- getMatrix(object)
  a <- matrix(0, nrow(m), ncol(m))
  a <- matrix(0, nrow(m), ncol(m))
  if (above) {
    a[m >= cutoff] <- 1
  } else {
    a[m <= cutoff] <- 1
  }
  diag(a) <- 0   # TODO decide if we should avoid self loops
  return(a)
}