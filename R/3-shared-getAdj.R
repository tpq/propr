#' Get Adjacency Matrix as indicated by permutation tests
#'
#' This function gets the significant pairs, according to the permutation tests. 
#' Then it fills the adjacency matrix with 1 if pair is significant, otherwise 0.
#' The diagonal is set to 0.
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate. Default is 0.05.
#' @param window_size An integer. Default is 1. When it is greater than 1,
#' the function will use a significant cutoff based on the moving
#' average of the FDR values. This is useful when the FDR values are
#' noisy and the user wants to smooth them out.
#' @param consider_negative_values A boolean or NULL. If NULL, default behaviour
#' is set for each metric. If TRUE, the function will consider both positive and
#' negative values. If FALSE, it will focus only on the positive values. This is 
#' relevant only for the \code{propr} object.
#' @return An adjacency matrix.
#'
#' @export
getAdjFDR <- function(object, fdr = 0.05, window_size = 1, consider_negative_values = NULL) {
  
  if (inherits(object, "propr")){
    adj <- getAdjFDR.propr(object, fdr=fdr, window_size=window_size, consider_negative_values=consider_negative_values)

  }else if(inherits(object, "propd")){
    adj <- getAdjFDR.propd(object, fdr=fdr, window_size=window_size)

  }else{
    stop("Please provide a 'propr' or 'propd' object.")
  }

  return(adj)
}

#' @rdname getAdjFDR
#' @section Methods:
#' \code{getAdjFDR.propd:}
#' This function creates an adjacency matrix with the significant
#' pairs from a \code{propd} object.
#' @export
getAdjFDR.propd <- 
  function(object, fdr = 0.05, window_size = 1) {

    # get cutoff
    cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size)

    # get theta matrix
    mat <- get_theta_matrix(object)

    # fill in significant pairs
    adj <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    rownames(adj) <- rownames(mat)
    colnames(adj) <- colnames(mat)
    adj[mat <= cutoff] <- 1

    return(adj)
  }

#' @rdname getAdjFDR
#' @section Methods:
#' \code{getAdjFDR.propr:}
#' This function creates an adjacency matrix with the significant
#' pairs from a \code{propd} object.
#' @export
getAdjFDR.propr <-
  function(object, fdr = 0.05, window_size = 1, consider_negative_values = NULL) {

    if (is.null(consider_negative_values)) consider_negative_values <- object@has_meaningful_negative_values
    
    if (!object@has_meaningful_negative_values & consider_negative_values) {
      message("Alert: negative values may not be relevant for this metric.")
    }
    if (object@has_meaningful_negative_values & !consider_negative_values) {
      message("Alert: try to set consider_negative_values to TRUE as the negative values may be relevant.")
    }

    # create empty matrix
    adj <- matrix(0, nrow = ncol(object@matrix), ncol = ncol(object@matrix))
    rownames(adj) <- rownames(object@matrix)
    colnames(adj) <- colnames(object@matrix)

    # for metrics with direct correlation
    if (object@direct){
      adj[
        object@matrix >= 0 &
        object@matrix >= getCutoffFDR(object, fdr=fdr, window_size=window_size)
      ] <- 1
      if (consider_negative_values) adj[
        object@matrix < 0 & 
        object@matrix <= getCutoffFDR(object, fdr=fdr, window_size=window_size, positive = FALSE)
      ] <- 1

    # for other metrics
    } else {
      adj[
        object@matrix >= 0 &
        object@matrix <= getCutoffFDR(object, fdr=fdr, window_size=window_size)
      ] <- 1
      if (consider_negative_values) adj[
        object@matrix < 0 & 
        object@matrix >= getCutoffFDR(object, fdr=fdr, window_size=window_size, positive = FALSE)
      ] <- 1
    }

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
getAdjFstat <- 
  function(object, pval = 0.05, fdr = TRUE) {
  
  if (!"Fstat" %in% colnames(object@results)) {
    stop("Please run updateF() on propd object before.")
  }

  if (pval < 0 | pval > 1) {
    stop("Provide a p-value cutoff from [0, 1].")
  }

  # get matrix
  mat <- get_theta_matrix(object)

  # create empty adjacency matrix
  adj <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  rownames(adj) <- rownames(mat)
  colnames(adj) <- colnames(mat)

  # fill in significant pairs
  cutoff <- getCutoffFstat(object, pval, fdr = fdr)
  adj[mat <= cutoff] <- 1

  return(adj)  
}
