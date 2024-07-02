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
getAdjFDR <- function(object, fdr = 0.05, positive = TRUE, negative = FALSE){
  results <- getSignificantResultsFDR(object, fdr = fdr, positive = positive, negative = negative)
  adj <- results2matRcpp(results, n = ncol(object@counts), diagonal = 0)
  return(adj)
}

getAdjFstat <- function(object, pval = 0.05, fdr = TRUE){
  results <- getSignificantResultsFstat(object, pval = pval, fdr = fdr)
  adj <- results2matRcpp(results, n = ncol(object@counts), diagonal = 0)
  return(adj)
}
