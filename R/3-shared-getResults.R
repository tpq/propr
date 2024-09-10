#' Get Results from Object
#'
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propr} or \code{propd} object.
#'
#' @param object A \code{propr} or \code{propd} object.
#'
#' @return A \code{data.frame} of results.
#'
#' @export
getResults <-
  function(object) {
    results <- object@results
    names <- colnames(object@counts)
    results$Partner <- names[results$Partner]
    results$Pair <- names[results$Pair]
    return(results)
  }

#' Get Significant Results from Object based on the permutation tests.
#'
#' This function retrieves results from a \code{propr} or \code{propd} object keeping only the
#' statistically significant pairs. The significance is determined by the cutoff value for which
#' the false discovery rate (FDR) is less or equal than the given value 'fdr'. The significant 
#' pairs are those that have a value greater/less or equal than the cutoff, depending on the metric.
#'
#' @inheritParams getAdjacencyFDR
#' @return A \code{data.frame} of results.
#'
#' @export
getSignificantResultsFDR <-
  function(object, fdr = 0.05, window_size = 1, tails = c('right', 'both')) {

    if (inherits(object, "propr")) {
      results <- getSignificantResultsFDR.propr(object, fdr=fdr, window_size=window_size, tails=tails)

    } else if(inherits(object, "propd")) {
      results <- getSignificantResultsFDR.propd(object, fdr=fdr, window_size=window_size)

    } else {
      stop("Please provide a 'propr' or 'propd' object.")
    }

    return(results)
}

#' @rdname getSignificantResultsFDR
#' @section Methods:
#' \code{getSignificantResultsFDR.propd:}
#' This function retrieves results from a \code{propd} object keeping 
#' only the statistically significant pairs.
#' @export
getSignificantResultsFDR.propd <- 
  function(object, fdr = 0.05, window_size = 1) {

    results <- getResults(object)
    cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size)
    if (cutoff) {
      return(results[which(results$theta <= cutoff), ])
    } else {
      return(results[0,])
    }
}

#' @rdname getSignificantResultsFDR
#' @section Methods:
#' \code{getSignificantResultsFDR.propr:}
#' This function retrieves results from a \code{propr} object keeping 
#' only the statistically significant pairs.
#' @export
getSignificantResultsFDR.propr <- 
  function(object, fdr = 0.05, window_size = 1, tails = c('right', 'both')) {

    # handle tails argument
    tails <- match.arg(tails)
    if (tails == 'both' & !object@direct) {
      warning("Running tails='right' instead")  
      tails <- 'right'   # set to right, when non negative values can be expected for a given metric.
    }

    # function to subset the results data frame based on the cutoff
    subsetBeyondCutoff <- function(data, cutoff) {
      if (!cutoff) return(data[0,])  # return empty data frame when no cutoff found
      if (object@direct) {
        return(data[which(abs(data$propr) >= abs(cutoff)), ])
      } else {
        return(data[which(data$propr <= cutoff), ])
      }
    }

    # correct fdr when considering both tails
    # so that the total fdr is the given value
    if (tails == 'both') fdr <- fdr / 2

    # define results data frame
    df <- getResults(object)
    results <- data.frame(matrix(ncol = ncol(df), nrow = 0))
    colnames(results) <- colnames(df)

    # get the significant positive values
    cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size)
    part <- df[which(df$propr >= 0),]
    results <- rbind(results, subsetBeyondCutoff(part, cutoff))

    # get the significant negative values
    if (tails == 'both') {
      cutoff <- getCutoffFDR(object, fdr=fdr, window_size=window_size, positive=F)
      part <- df[which(df$propr < 0),]
      results <- rbind(results, subsetBeyondCutoff(part, cutoff))
    }

    return(results)
}

#' Get Significant Results based on the F-stats.
#'
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propd} object keeping only the statistically 
#'  significant pairs. Note that it can only be applied to theta_d,
#'  as updateF only works for theta_d.
#'
#' @inheritParams getAdjacencyFstat
#' @return A \code{data.frame} of results.
#'
#' @export
getSignificantResultsFstat <- 
  function(object, pval = 0.05, fdr_adjusted = TRUE) {

    if (!"Fstat" %in% colnames(object@results)) {
      stop("Please run updateF() on propd object before.")
    }
    if (pval < 0 | pval > 1) {
      stop("Provide a p-value cutoff from [0, 1].")
    }

    # get results data frame
    results <- getResults(object)

    # get significant theta based on the FDR adjusted empirical p-values
    if (fdr_adjusted) {
      message("Alert: Returning the significant pairs based on the FDR adjusted p-values.")
      results <- results[which(results$FDR <= pval), ]

    # get siniginicant theta based on the F-statistic cutoff
    } else {
      message("Alert: Returning the significant pairs based on the F-statistic cutoff.")
      cutoff <- getCutoffFstat(object, pval = pval, fdr_adjusted = FALSE)
      if (cutoff) {
        results <- results[which(results$theta <= cutoff), ]
      } else {
        results <- results[0,]
      }
    }

    return(results)
}
