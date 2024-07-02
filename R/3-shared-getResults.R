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
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propr} or \code{propd} object keeping only the
#'  statistically significant pairs.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param fdr A numeric. The false discovery rate to use.
#' @param positive A boolean. In the case of \code{propr} object,
#' toggles whether to return the significant positive values (>=0).
#' @param negative A boolean. In the case of \code{propr} object,
#' toggles whether to return the significant negative values (<0).
#' @return A \code{data.frame} of results.
#'
#' @export
getSignificantResultsFDR <-
  function(object, fdr = 0.05, positive = TRUE, negative = FALSE) {

    if(inherits(object, "propr")){
      if (!positive & !negative) stop("Please provide either positive or negative.")
      results <- getSignificantResultsFDR.propr(object, fdr = fdr, positive = positive, negative = negative)
    }else if(inherits(object, "propd")){
      results <- getSignificantResultsFDR.propd(object, fdr = fdr)
    }else{
      stop("Please provide a 'propr' or 'propd' object.")
    }
    return(results)
}

getSignificantResultsFDR.propd <- 
  function(object, fdr = 0.05) {
    results <- getResults(object)
    cutoff <- getCutoffFDR(object, fdr = fdr)
    results <- results[which(results$theta <= cutoff), ]
    return(results)
}

getSignificantResultsFDR.propr <- 
  function(object, fdr = 0.05, positive = TRUE, negative = FALSE) {

    # get metric
    metric <- object@metric

    # define results data frame
    df <- getResults(object)
    results <- data.frame(matrix(ncol = ncol(df), nrow = 0))
    colnames(results) <- colnames(df)

    # get the significant positive values
    if (positive) {
      cutoff <- getCutoffFDR(object, fdr = fdr, positive = TRUE)
      part <- df[which(df$propr >= 0),]
      if (metric_is_direct(metric)) {
        tmp <- part[which(part$propr >= cutoff), ]
      } else {
        tmp <- part[which(part$propr <= cutoff), ]
      }
      results <- rbind(results, tmp)
    }

    # get the significant negative values
    if (negative) {
      cutoff <- getCutoffFDR(object, fdr = fdr, positive = FALSE)
      part <- df[which(df$propr < 0),]
      if (metric_is_direct(metric)) {
        tmp <- df[which(df$propr <= cutoff), ]
      } else {
        tmp <- df[which(df$propr >= cutoff), ]
      }
      results <- rbind(results, tmp)
    }

    return(results)
}

#' Get Significant Results based on the F-stats.
#'
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propd} object keeping only the statistically 
#'  significant pairs. Note that it can only be applied to theta_d.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param pval A numeric. The p-value to use.
#' @param fdr A boolean. If TRUE, the function will return the
#' significant pairs based on the FDR adjusted p-values. Otherwise,
#' it will return the significant pairs based on the F-statistic cutoff.
#' @return A \code{data.frame} of results.
#'
#' @export
getSignificantResultsFstat <- 
  function(object, pval = 0.05, fdr = TRUE) {

    if (!"Fstat" %in% colnames(object@results)) {
      stop("Please run updateF() on propd object before.")
    }

    if (pval < 0 | pval > 1) {
      stop("Provide a p-value cutoff from [0, 1].")
    }

    # get results data frame
    results <- getResults(object)

    # get significant theta based on the FDR adjusted empirical p-values
    if (fdr) {
      message("Alert: Returning the significant pairs based on the FDR adjusted p-values.")
      results <- results[which(results$FDR <= pval), ]

    # get siniginicant theta based on the F-statistic cutoff
    }else{
      message("Alert: Returning the significant pairs based on the F-statistic cutoff.")
      cutoff <- getCutoffFstat(object, pval = pval, fdr = FALSE)
      results <- results[which(results$Fstat <= cutoff), ]
    }

    return(results)
}
