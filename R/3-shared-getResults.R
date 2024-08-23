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
#' @param consider_negative_values A boolean. If TRUE, the function will
#' consider also the negative values. Otherwise, it will focus only on the
#' positive values. This is relevant only for the \code{propr} object.
#' @return A \code{data.frame} of results.
#'
#' @export
getSignificantResultsFDR <-
  function(object, fdr = 0.05, consider_negative_values = FALSE) {

    if(inherits(object, "propr")){
      results <- getSignificantResultsFDR.propr(object, fdr, consider_negative_values)

    }else if(inherits(object, "propd")){
      results <- getSignificantResultsFDR.propd(object, fdr)
    }else{

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
  function(object, fdr = 0.05) {

    results <- getResults(object)
    cutoff <- getCutoffFDR(object, fdr = fdr)
    results <- results[which(results$theta <= cutoff), ]

    return(results)
}

#' @rdname getSignificantResultsFDR
#' @section Methods:
#' \code{getSignificantResultsFDR.propr:}
#' This function retrieves results from a \code{propr} object keeping 
#' only the statistically significant pairs.
#' @export
getSignificantResultsFDR.propr <- 
  function(object, fdr = 0.05, consider_negative_values = FALSE) {

    if (!object@has_meaningful_negative_values & consider_negative_values) {
      message("Alert: negative values may not be relevant for this metric.")
    }
    if (object@has_meaningful_negative_values & !consider_negative_values) {
      message("Alert: try to set consider_negative_values to TRUE as the negative values may be relevant.")
    }

    # function to subset the results data frame based on the cutoff
    subsetBeyondCutoff <- function(data, cutoff) {
      if (object@direct) {
        data <- data[which(abs(data$propr) >= abs(cutoff)), ]
      } else {
        data <- data[which(abs(data$propr) <= abs(cutoff)), ]
      }
      return(data)
    }

    # define results data frame
    df <- getResults(object)
    results <- data.frame(matrix(ncol = ncol(df), nrow = 0))
    colnames(results) <- colnames(df)

    # get the significant positive values
    cutoff <- getCutoffFDR(object, fdr = fdr, positive = TRUE)
    part <- df[which(df$propr >= 0),]
    results <- rbind(results, subsetBeyondCutoff(part, cutoff))

    # get the significant negative values
    if (consider_negative_values) {
      cutoff <- getCutoffFDR(object, fdr = fdr, positive = FALSE)
      part <- df[which(df$propr < 0),]
      results <- rbind(results, subsetBeyondCutoff(part, cutoff))
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
