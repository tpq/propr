#' Get Adjacency Matrix as indicated by permutation tests
#'
#' This function gets the significant pairs, according to the permutation tests. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' The diagonal is set to 0.
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate. Default is 0.05.
#' @param positive A boolean. In the case of \code{propr} object,
#' toggles where to retain the significant pairs with positive values (>=0).
#' @param negative A boolean. In the case of \code{propr} object,
#' toggles where to retain the significant pairs with negative values (<0).
#'
#' @return An adjacency matrix.
#'
#' @export
getAdjFDR <- function(object, fdr = 0.05, positive = TRUE, negative = FALSE){
  results <- getSignificantResultsFDR(object, fdr = fdr, positive = positive, negative = negative)
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
