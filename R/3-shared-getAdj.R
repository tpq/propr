#' Get Adjacency Matrix as indicated by permutation tests
#'
#' This function gets the significant pairs, according to the permutation tests. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' The diagonal is set to 0.
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate. Default is 0.05.
#'
#' @return An adjacency matrix.
#'
#' @export
getAdjFDR <- function(object, fdr = 0.05){
  results <- getSignificantResultsFDR(object, fdr = fdr)
  adj <- results2matRcpp(results, n = ncol(object@counts), diagonal = 0)
  return(adj)
}

#' Get Adjacency Matrix as indicated by F-statistics
#'
#' This function gets the significant pairs, according to the F-statistics. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' The diagonal is set to 0.
#' @param object A \code{propd} or \code{propr} object.
#' @param pval A float value for the p-value. Default is 0.05.
#' @param fdr A boolean. If TRUE, use the the FDR- adjusted p-values. Otherwise,
#' get significant pairs based on the theoretical F-statistic cutoff.
#' @return An adjacency matrix.
#'
#' @export
getAdjFstat <- function(object, pval = 0.05, fdr = TRUE){
  results <- getSignificantResultsFstat(object, pval = pval, fdr = fdr)
  adj <- results2matRcpp(results, n = ncol(object@counts), diagonal = 0)
  return(adj)
}
