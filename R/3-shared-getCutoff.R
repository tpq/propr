#' Get a meaningful cutoff based on the FDR values from permutation tests.
#' 
#' @param object A \code{propd} or \code{propr} object.
#' @param fdr A float value for the false discovery rate.
#' Default is 0.05.
#' @param window_size An integer. Default is 1. When it is greater than 1,
#' the function will return a meaningful cutoff based on the moving
#' average of the FDR values. This is useful when the FDR values are
#' noisy and the user wants to smooth them out.
#' @return A cutoff value.
#' @export
getCutoffFDR <- function(object, fdr = 0.05, window_size = 1) {
  NVTX_PUSH("getCutoffFDR", 0)

  if (!"fdr" %in% slotNames(object)) {
    NVTX_POP()
    stop("Please run updateCutoffs() before calling this function.")
  }
  if (nrow(object@fdr) == 0) {
    NVTX_POP()
    stop("No FDR values found. Please run updateCutoffs() before calling this function.")
  }
  if (fdr < 0 | fdr > 1) {
    NVTX_POP()
    stop("Provide a FDR cutoff from [0, 1].")
  }

  # get data frame
  df <- object@fdr

  # apply moving average to FDR values, if window_size > 1
  if (window_size > 1) {
    NVTX_PUSH("moving_average_FDR", 0)
    message("Applying moving average to FDR values.")
    df$FDR <- getMovingAverage(df$FDR, window_size)
    NVTX_POP()
  }

  # get index of FDR values below the threshold
  NVTX_PUSH("threshold_index_FDR", 0)
  index <- (df$FDR <= fdr) & (is.finite(df$FDR))
  NVTX_POP()
  if (!any(index)) {
    warning("No significant cutoff found for the given FDR = ", fdr)
    NVTX_POP()  # getCutoffFDR
    return(FALSE)
  }

  # get cutoff
  NVTX_PUSH("compute_cutoff_FDR", 0)
  direct <- FALSE
  if (inherits(object, "propr")) direct <- object@direct
  if (direct) {
    cutoff <- min(df$cutoff[index])
  } else {
    cutoff <- max(df$cutoff[index])
  }
  NVTX_POP()  # compute_cutoff_FDR

  NVTX_POP()  # getCutoffFDR
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
#' @param fdr_adjusted A boolean. Toggles whether to calculate the theta
#' cutoff for an FDR-adjusted p-value.
#' @return A cutoff of theta from [0, 1].
#' @export
getCutoffFstat <- function(object, pval = 0.05, fdr_adjusted = FALSE) {
  NVTX_PUSH("getCutoffFstat", 0)

  if (!"Fstat" %in% colnames(object@results)) {
    NVTX_POP()
    stop("Please run updateF() on propd object before.")
  }
  if (pval < 0 | pval > 1) {
    NVTX_POP()
    stop("Provide a p-value cutoff from [0, 1].")
  }

  if (fdr_adjusted) {
    NVTX_PUSH("empiric_FDR_cutoff", 0)
    message("Alert: Returning an empiric cutoff based on the $FDR slot.")
    index <- (object@results$FDR <= pval) & (is.finite(object@results$FDR))
    if (any(index)) {
      cutoff <- max(object@results$theta[index])
    } else {
      warning("No significant cutoff found for the given p-value.")
      cutoff <- FALSE
    }
    NVTX_POP()  # empiric_FDR_cutoff
  } else {
    NVTX_PUSH("theoretical_Fstat_cutoff", 0)
    message("Alert: Returning an cutoff based on the F-statistic.")
    # Compute based on theory
    K <- length(unique(object@group))
    N <- length(object@group) + object@dfz # population-level metric (i.e., N)
    Q <- stats::qf(pval, K - 1, N - K, lower.tail = FALSE)
    # theta = (N-2)/(Q+(N-2))
    cutoff <- (N - 2) / (Q + (N - 2))
    NVTX_POP()  # theoretical_Fstat_cutoff
  }

  NVTX_POP()  # getCutoffFstat
  return(cutoff)
}


#' Caclulate the moving average of a vector.
#' @param values A numeric vector.
#' @param window_size An integer. The size of the window to calculate the
#' moving average. Default is 1.
getMovingAverage <- function(values, window_size = 1) {
  NVTX_PUSH("getMovingAverage", 0)

  if (any(is.na(values))) {
    message("Moving averages are calculated for a vector containing NAs.")
  }

  # Initialize the result vector
  NVTX_PUSH("init_result", 0)
  n <- length(values)
  result <- numeric(n)
  NVTX_POP()  # init_result
  
  NVTX_PUSH("moving_average_loop", 0)
  for (i in 1:n) {
    # Determine the window indices
    if (window_size %% 2 == 0) {
      start_idx <- max(1, i - (window_size / 2 - 1))
      end_idx <- min(n, i + window_size / 2)
    } else {
      start_idx <- max(1, i - floor(window_size / 2))
      end_idx <- min(n, i + floor(window_size / 2))
    }
    
    # Calculate the average for the current window
    if (is.finite(values[i])) {
      # NA values are removed, to avoid propagation of NAs
      result[i] <- mean(values[start_idx:end_idx], na.rm = TRUE)
    } else {
      # this keeps the NA (or non-finite) values corresponding to that position
      result[i] <- values[i]
    }
  }
  NVTX_POP()  # moving_average_loop
  
  NVTX_POP()  # getMovingAverage
  return(result)
}
