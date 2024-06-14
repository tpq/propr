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

  # get gene-gene matrix
  m <- getMatrix(object)

  # create adjacency matrix
  a <- matrix(0, nrow(m), ncol(m))
  rownames(a) <- rownames(m)
  colnames(a) <- colnames(m)

  # fill in values
  if (above) {
    a[m >= cutoff] <- 1
  } else {
    a[m <= cutoff] <- 1
  }
  diag(a) <- 0
  
  return(a)
}