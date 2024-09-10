#' Get Adjacency Matrix as indicated by permutation tests.
#' 
#' This function gets the significant pairs according to the permutation tests. Then it fills 
#' the adjacency matrix with 1 if pair is significant, otherwise 0. The significance is determined
#' by the cutoff value for which the false discovery rate (FDR) is less or equal than the given 
#' value 'fdr'. The significant pairs are those that have a value greater/less or equal than the 
#' cutoff, depending on the metric.
#'
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate. Default is 0.05.
#' @param window_size An integer. Default is 1. When it is greater than 1, the FDR
#' values would be smoothed out by a moving average of the given window size.
#' @param tails 'right' or 'both'. 'right' is for one-sided on the right. 'both' is to 
#' combine one-sided on the right (positive values) and left (negative values). This is only 
#' relevant for \code{propr} objects, as \code{propd} objects are always one-sided and only 
#' have positive values.
#' @return An adjacency matrix.
#'
#' @export
getAdjacencyFDR <- 
  function(object, fdr = 0.05, window_size = 1, tails = c('right', 'both')) {

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
  function(object, fdr = 0.05, window_size = 1, tails = c('right', 'both')) {

    # handle tails argument
    tails <- match.arg(tails)
    if (tails == 'both' & !object@direct) {
      warning("Running tails='right' instead")  
      tails <- 'right'   # set to right, when non negative values can be expected for a given metric.
    }

    # correct fdr when considering both tails
    # so that the total fdr is the given value
    if (tails == 'both') fdr <- fdr / 2

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
        adj[object@matrix <= cutoff] <- 1  # don't have negative values
      }
    }

    # for negative tail
    if (tails == 'both'){
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
