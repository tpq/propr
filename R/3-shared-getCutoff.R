#' Get a meaningful cutoff based on the FDR values from permutation tests.
#' 
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate.
#' Default is 0.05.
#' @param positive A boolean. Toggles whether to return the cutoff with
#' positive value. Otherwise, return cutoff with negative value.
#' @return A cutoff value.
#' @export
getCutoffFDR <- function(object, fdr = 0.05, positive = TRUE) {
  if (!"fdr" %in% slotNames(object)) {
    stop("Please run updateCutoffs() before calling this function.")
  }
  
  if (fdr < 0 | fdr > 1) {
    stop("Provide a FDR cutoff from [0, 1].")
  }

  # subset data frame depending if considering positive or negative values
  df <- object@fdr
  if (positive) {
    df <- df[df$cutoff >= 0, ]
  } else {
    df <- df[df$cutoff < 0, ]
  }

  # get index of FDR values below the threshold
  index <- (df$FDR <= fdr) & (is.finite(df$FDR))
  if (!any(index)) {
    stop("No pairs below FDR.")  # TODO should we throw an error, or just return FALSE?
  }

  # get cutoff
  if ( (class(object) == 'propr') && (metric_is_direct(object@metric)) ) {
    cutoff <- min(df$cutoff[index])
  } else{
    cutoff <- max(df$cutoff[index])
  }

  return(cutoff)
}

#' Calculate a theta Cutoff based on the F-statistic.
#'
#' This function uses the F distribution to calculate a cutoff of
#'  theta for a p-value given by the \code{pval} argument.
#'
#' If the argument \code{fdr = TRUE}, this function returns the
#'  empiric cutoff that corresponds to the FDR-adjusted p-value
#'  stored in the \code{@@results$FDR} slot.
#'
#' @param object A \code{\link{propd}} object.
#' @param pval A p-value at which to calculate a theta cutoff.
#' @param fdr A boolean. Toggles whether to calculate the theta
#'  cutoff for an FDR-adjusted p-value.
#' @return A cutoff of theta from [0, 1].
#' @export
getCutoffFstat <- function(object, pval = 0.05, fdr = FALSE) {
  if (!"Fstat" %in% colnames(object@results)) {
    stop("Please run updateF() on propd object before calling qtheta.")
  }

  if (pval < 0 | pval > 1) {
    stop("Provide a p-value cutoff from [0, 1].")
  }

  if (fdr) {
    message("Alert: Returning an empiric cutoff based on the $FDR slot.")
    index <- (object@results$FDR < pval) & (is.finite(object@results$FDR))
    if (any(index)) {
      cutoff <- max(object@results$theta[index])
    } else{
      stop("No pairs below p-value.")
    }

  } else{
    message("Alert: Returning an cutoff based on the F-statistic.")
    # Compute based on theory
    K <- length(unique(object@group))
    N <-
      length(object@group) + object@dfz # population-level metric (i.e., N)
    Q <- stats::qf(pval, K - 1, N - K, lower.tail = FALSE)
    # # Fstat <- (N - 2) * (1 - object@theta$theta) / object@theta$theta
    # # Q = Fstat
    # # Q = (N-2) * (1-theta) / theta
    # # Q / (N-2) = (1/theta) - 1
    # # 1/theta = Q / (N-2) + 1 = Q(N-2)/(N-2)
    # # theta = (N-2)/(Q+(N-2))
    cutoff <- (N - 2) / (Q + (N - 2))
  }

  return(cutoff)
}
