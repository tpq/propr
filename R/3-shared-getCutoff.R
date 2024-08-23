#' Get a meaningful cutoff based on the FDR values from permutation tests.
#' 
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate.
#' Default is 0.05.
#' @param window_size An integer. Default is 1. When it is greater than 1,
#' the function will return a meaningful cutoff based on the moving
#' average of the FDR values. This is useful when the FDR values are
#' noisy and the user wants to smooth them out.
#' @param positive A boolean. In the case of \code{propr} object, 
#' toggles whether to return the cutoff for positive values (>=0). 
#' Otherwise, return cutoff for negative values (<0). Default is TRUE.
#' This argument is ignored for \code{propd} objects.
#' @return A cutoff value.
#' @export
getCutoffFDR <- function(object, fdr = 0.05, window_size = 1, positive = TRUE) {
  if (!"fdr" %in% slotNames(object)) {
    stop("Please run updateCutoffs() before calling this function.")
  }
  if (nrow(object@fdr) == 0) {
    stop("No FDR values found. Please run updateCutoffs() before calling this function.")
  }

  if (fdr < 0 | fdr > 1) {
    stop("Provide a FDR cutoff from [0, 1].")
  }

  if (inherits(object, "propr")) {
    return(getCutoffFDR.propr(object, fdr, window_size, positive))

  } else if (inherits(object, "propd")) {
    return(getCutoffFDR.propd(object, fdr, window_size))
    
  } else{
    stop("Provided 'object' not recognized.")
  }
}

#' @rdname getCutoffFDR
#' Same as \code{getCutoffFDR}, but for \code{propr} objects.
#' @export
getCutoffFDR.propr <- function(object, fdr = 0.05, window_size = 1, positive = TRUE) {

  # subset FDR data frame, to only include positive cutoffs
  df <- object@fdr
  if (positive) { df <- df[df$cutoff >= 0,] } else { df <- df[df$cutoff < 0,] }

  # apply moving average to FDR values
  if (window_size > 1) {
    message("Applying moving average to FDR values.")
    df$FDR <- getMovingAverage(df$FDR, window_size)
  }

  # get index of FDR values below the threshold
  index <- (df$FDR <= fdr) & (is.finite(df$FDR))
  if (!any(index)) {
    stop("No pairs below FDR.")  # TODO should we throw an error, or just return FALSE?
  }

  # get cutoff
  if (object@direct) {
    cutoff <- min(abs(df$cutoff[index]))
  } else{
    cutoff <- max(abs(df$cutoff[index]))
  }
  if (!positive) cutoff <- -cutoff

  return(cutoff)
}

#' @rdname getCutoffFDR
#' Same as \code{getCutoffFDR}, but for \code{propd} objects.
#' @export
getCutoffFDR.propd <- function(object, fdr = 0.05, window_size = 1) {

  # get df
  df <- object@fdr
  if (any(df$cutoff < 0)) stop("Negative cutoffs were found for theta values. Please check")

  # apply moving average to FDR values
  if (window_size > 1) {
    message("Applying moving average to FDR values.")
    df$FDR <- getMovingAverage(df$FDR, window_size)
  }
  
  # get index of FDR values below the threshold
  index <- (df$FDR <= fdr) & (is.finite(df$FDR))
  if (!any(index)) {
    stop("No pairs below FDR.")  # TODO should we throw an error, or just return FALSE?
  }

  # get cutoff
  cutoff <- max(df$cutoff[index])
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
#' cutoff for an FDR-adjusted p-value.
#' @return A cutoff of theta from [0, 1].
#' @export
getCutoffFstat <- function(object, pval = 0.05, fdr = FALSE) {
  if (!"Fstat" %in% colnames(object@results)) {
    stop("Please run updateF() on propd object before.")
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
    N <- length(object@group) + object@dfz # population-level metric (i.e., N)
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

#' Caclulate the moving average of a vector.
#' @param values A numeric vector.
#' @param window_size An integer. The size of the window to calculate the
#' moving average. Default is 1.
getMovingAverage <- function(values, window_size = 1) {

  if (any(is.na(values))) {
    message("Moving averages are calculated for a vector containing NAs.")
  }

  # Initialize the result vector
  n <- length(values)
  result <- numeric(n)
  
  for (i in 1:n) {
    # Determine the window indices
    if (window_size %% 2 == 0) {
      start_idx <- max(1, i - (window_size / 2 - 1))
      end_idx <- min(n, i + window_size / 2)
    }else{
      start_idx <- max(1, i - floor(window_size / 2))
      end_idx <- min(n, i + floor(window_size / 2))
    }
    
    # Calculate the average for the current window
    if (is.finite(values[i])){
      result[i] <- mean(values[start_idx:end_idx], na.rm=TRUE)  # NA values are removed, to avoid propagation of NAs
    }else{
      result[i] <- values[i]  # this keeps the NA values corresponding to that position
    }
  }
  
  return(result)
}