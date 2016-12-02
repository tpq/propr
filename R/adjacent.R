#' Make Adjacency Object
#'
#' This function uses pairs indexed in the \code{@@pairs}
#'  slot to build a symmetric adjacency matrix.
#'
#' @inheritParams bucket
#'
#' @return A \code{propr} object with the adjacency matrix
#'  saved to the \code{@@matrix} slot.
#'
#' @export
adjacent <- function(rho){

  if(length(rho@pairs) == 0 | class(rho) != "propr"){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  N <- nrow(rho@matrix)
  mat <- matrix(0, N, N)
  mat[rho@pairs] <- 1
  diag(mat) <- 1
  symRcpp(mat)

  adj <- rho
  adj@matrix <- mat

  return(adj)
}
