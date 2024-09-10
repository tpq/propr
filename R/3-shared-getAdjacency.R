#' Get Adjacency Matrix as indicated by permutation tests
#'
#' This function gets the significant pairs, according to the permutation tests. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate. Default is 0.05.
#' @param window_size An integer. Default is 1. When it is greater than 1, the FDR
#' values would be smoothed out by a moving average of the given window size.
#' @param tails 1 for one-sided on the right, and 2 for two-sided. When NULL, use default
#' test for the given metric. This is only relevant for \code{propr} objects.
#' @return An adjacency matrix.
#'
#' @export
getAdjacencyFDR <- function(object, fdr = 0.05, window_size = 1, tails = NULL) {
  
  if (inherits(object, "propr")){
    adj <- getAdjacencyFDR.propr(object, fdr=fdr, window_size=window_size, tails=tails)

  }else if(inherits(object, "propd")){
    adj <- getAdjacencyFDR.propd(object, fdr=fdr, window_size=window_size)

  }else{
    stop("Please provide a 'propr' or 'propd' object.")
  }

  return(adj)
}

#' @rdname getAdjacencyFDR
#' @section Methods:
#' \code{getAdjacencyFDR.propd:}
#' This function creates an adjacency matrix with the significant
#' pairs from a \code{propd} object.
#' @export
getAdjacencyFDR.propd <- 
  function(object, fdr = 0.05, window_size = 1) {

    # get cutoff
    cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size)

    # get theta matrix
    mat <- getMatrix(object)

    # create empty matrix
    adj <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    rownames(adj) <- rownames(mat)
    colnames(adj) <- colnames(mat)
    diag(adj) <- 1

    # fill in significant pairs
    if (cutoff) adj[mat <= cutoff] <- 1

    return(adj)
  }

#' @rdname getAdjacencyFDR
#' @section Methods:
#' \code{getAdjacencyFDR.propr:}
#' This function creates an adjacency matrix with the significant
#' pairs from a \code{propd} object.
#' @export
getAdjacencyFDR.propr <-
  function(object, fdr = 0.05, window_size = 1, tails = NULL) {

    if (is.null(tails)) {
      tails <- ifelse(object@has_meaningful_negative_values, 2, 1)
    }
    if (tails != 1 & tails != 2) {
        stop("Please provide a valid value for tails: 1 or 2.")
    }
    if (tails == 1 & object@has_meaningful_negative_values) {
      warning("Significant pairs are chosen based on one-sided FDR test.")
    }
    if (tails == 2 & !object@direct) {
      stop("Two-sided FDR is not available for this metric.")
    }

    # correct fdr when two-sided
    if (tails == 2) fdr <- fdr / 2

    # create empty matrix
    adj <- matrix(0, nrow = ncol(object@matrix), ncol = ncol(object@matrix))
    rownames(adj) <- rownames(object@matrix)
    colnames(adj) <- colnames(object@matrix)
    diag(adj) <- 1

    # for positive tail
    cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size)
    if (cutoff) {
      if (object@direct) {
        adj[object@matrix >= 0 & object@matrix >= cutoff] <- 1
      } else {
        adj[object@matrix <= cutoff] <- 1
      }
    }

    # for negative tail
    if (tails == 2){
      cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size, positive=F)
      if (cutoff) adj[object@matrix < 0 & object@matrix <= cutoff] <- 1
    }

    return(adj)
  }

#' Get Adjacency Matrix as indicated by F-statistics
#'
#' This function gets the significant pairs, according to the F-statistics. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' Note that it can only be applied to theta_d, as updateF only works for theta_d.
#' 
#' @param object A \code{propd} or \code{propr} object.
#' @param pval A float value for the p-value. Default is 0.05.
#' @param fdr_adjusted A boolean. If TRUE, use the the FDR- adjusted p-values. 
#' Otherwise, get significant pairs based on the theoretical F-statistic cutoff.
#' @return An adjacency matrix.
#'
#' @export
getAdjacencyFstat <- 
  function(object, pval = 0.05, fdr_adjusted = TRUE) {
  
  # get matrix
  mat <- getMatrix(object)

  # create empty adjacency matrix
  adj <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  rownames(adj) <- rownames(mat)
  colnames(adj) <- colnames(mat)

  # fill in significant pairs
  cutoff <- getCutoffFstat(object, pval, fdr_adjusted = fdr_adjusted)
  if (cutoff) adj[mat <= cutoff] <- 1

  return(adj)  
}
