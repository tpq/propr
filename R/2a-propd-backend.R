#' Calculate Theta and Related Statistics
#'
#' This function calculates theta and related statistics based on the input
#'  count matrix and other parameters. The function provides various options
#'  for calculating theta (theta_d, theta_e, theta_f, theta_g).
#'
#' @inheritParams propd
#' @param lrv If LRV is provided, it is not computed within the function.
#' @param only A character vector specifying the type of theta to calculate.
#' @param weighted A logical value indicating whether weighted calculations
#'  should be performed. 
#' @param shrink A logical value indicating whether to apply shrinkage
#'
#' @return A data frame containing the computed theta values and
#'  related statistics, depending on the `only` parameter.
#'
#' @examples
#' # Sample input count data and group assignments
#' data <- iris[1:100, 1:4]
#' group <- iris[1:100, 5]
#'
#' # Calculate all theta types
#' result_all <- calculate_theta(data, group, alpha = 0.5)
#'
#' # Calculate only theta_d
#' result_theta_d <- calculate_theta(data, group, alpha = 0.5, only = "theta_d")
#'
#' @export
calculate_theta <-
  function(counts,
           group,
           alpha = NA,
           lrv = NA,
           only = "all",
           weighted = FALSE,
           shrink = FALSE) {

    # count matrix
    ct <- as.matrix(counts)

    # define recompute lrv or not
    if (identical(lrv, NA)) {
      firstpass <- TRUE
    } else{
      firstpass <- FALSE
    }

    # Get groups
    groups <- lapply(unique(group), function(g)
      g == group)
    ngrp <- length(unique(group))

    ##############################################################################
    ### Use weights
    ##############################################################################

    # calculate weights using limma sample weights
    if (weighted) {
      message("Alert: Calculating Limma's reliability weights for samples.")
      packageCheck("limma")
      
      #use clr-transform of the counts for quality weights from limma: 
      design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))

      # limma requires a pseudocount of +1 if zeros are present
      if (any(counts == 0)) {    
        X <- as.matrix(counts + 1)  
      } else{
        X <- as.matrix(counts)
      }

      # calculate geometric mean
      logX <- log(X)
      z.geo <- rowMeans(logX)

      # calculate clr-transformed data
      z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
      # scale counts by mean of geometric mean, this rescales data to similar magnitude as original counts
      lz.sr <- t(z.lr + mean(z.geo)) #corresponds to log(z.sr) in updateF function

      # use quality weights from limma:
      aw <- limma::arrayWeights(lz.sr, design) 
      weights <- t(sweep(matrix(1, nrow(lz.sr), ncol(lz.sr)), 2, aw, `*`)) #get the correct dimensions

      if (nrow(weights) != nrow(counts) | ncol(weights) != ncol(counts)) {
        stop("The matrix dimensions of 'weights' must match the matrix dimensions 'counts'.")
      }
    }

    # use weights for lrv modifiers, if provided
    if (weighted) {
      W <- weights
      ps <- lapply(groups, function(g) omega(W[g,]))
      names(ps) <- paste0("p", 1:ngrp)
      p <- omega(W)
    } else {
      W <- ct
      ps <- lapply(groups, function(g) sum(g) - 1)
      names(ps) <- paste0("p", 1:ngrp)
      p <- length(group) - 1
    }

    ##############################################################################
    ### Handle zeros
    ##############################################################################

    # Replace zeros if any
    # Logratio and theta cannot be computed with zeros
    if (any(ct == 0) && !is.na(alpha)) {
      ct <- simple_zero_replacement(ct)
    }

    ##############################################################################
    ### Calculate logratio variance and differential proportionality theta
    ##############################################################################

    # calculate logratio variance based on shrunk covariance matrix
    if (shrink) {
      if (weighted) {
        stop("Shrinkage is not available for weighted computation yet.")
      }
      if (firstpass)
        lrv <- lrv_with_shrinkage(ct)
      lrvs <- 
        lapply(groups, function(g)
          lrv_with_shrinkage(ct[g,]))
      names(lrvs) <- paste0("lrv", 1:ngrp)
    } else {

      # Calculate weighted and/or alpha-transformed LRVs -- W not used if weighted = FALSE
      if (firstpass)
        lrv <- lrv(ct, W, weighted, alpha, ct, W)
      lrvs <-
        lapply(groups, function(g)
          lrv(ct[g,], W[g,], weighted, alpha, ct, W))
      names(lrvs) <- paste0("lrv", 1:ngrp)
    }

    # Calculate LRM (using alpha-based LRM if appropriate)
    if (only == "all") {
      if (firstpass)
        lrm <- lrm(ct, W, weighted, alpha, ct, W)
      lrms <-
        lapply(groups, function(g)
          lrm(ct[g,], W[g,], weighted, alpha, ct, W))
      names(lrms) <- paste0("lrm", 1:ngrp)
    }

    # Replace NaN thetas (from VLR = 0 or VLR = NaN) with 1
    lrv0 <-
      Reduce("|", lapply(lrvs, is.na)) |
      is.na(lrv) | (lrv == 0) # aVLR triggers NaN
    replaceNaNs <- any(lrv0)
    if (replaceNaNs) {
      if (firstpass)
        message("Alert: Replacing NaN theta values with 1.")
    }

    # Calculate within-group sums-of-squares (used to calculate theta)
    SS <- lapply(1:ngrp, function(i)
      ps[[i]] * lrvs[[i]])

    ##############################################################################
    ### Build data.frame of results with computed theta
    ##############################################################################

    # Build all theta types unless only != "all"
    if (only == "all" | only == "theta_d") {
      theta <- Reduce(`+`, SS) / (p * lrv)
      if (replaceNaNs)
        theta[lrv0] <- 1
      if (only == "theta_d")
        return(theta)
    }

    if (only == "all" | only == "theta_e") {
      theta_e <- 1 - Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_e[lrv0] <- 1
      if (only == "theta_e")
        return(theta_e)
    }

    if (only == "all" | only == "theta_f") {
      theta_f <- Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_f[lrv0] <- 1
      if (only == "theta_f")
        return(theta_f)
    }

    if (only == "all" | only == "theta_g") {
      theta_g <- Reduce("pmin", SS) / (p * lrv)
      if (replaceNaNs)
        theta_g[lrv0] <- 1
      if (only == "theta_g")
        return(theta_g)
    }

    labels <- labRcpp(ncol(counts))
    return(
      data.frame(
        "Partner" = labels[[1]],
        "Pair" = labels[[2]],
        "theta" = theta,
        "theta_e" = theta_e,
        "theta_f" = theta_f,
        "theta_g" = theta_g,
        "lrv" = lrv,
        lrvs,
        "lrm" = lrm,
        lrms,
        "p" = p,
        ps
      )
    )
  }

#' Calculate Logratio Variance with shrinkage
#'
#' This function computes the logratio variance (LRV) with the option
#'  to apply shrinkage. It uses the `corpcor` package to compute a shrunk
#'  covariance matrix and then converts it to a logratio variance matrix.
#'
#' @param ct A count matrix.
#' @param shrink A logical value indicating whether to apply shrinkage.
#' @return A shrunk logratio variance matrix.
lrv_with_shrinkage <- function(ct, shrink = TRUE) {
  
  # compute covariance matrix on the log data
  if (shrink) {
    cov_matrix <- corpcor::cov.shrink(log(ct))
  } else {
    cov_matrix <- stats::cov(log(ct))
  }

  # convert shrunked covariance matrix to shrunked logratio variance matrix
  diag_matrix <- diag(cov_matrix)
  outer_sum <- outer(diag_matrix, diag_matrix, "+")
  lrv <- outer_sum - 2 * cov_matrix

  # it is symmetric
  # we get the upper triangle, so that it is coherent with the function propr:::lrv
  lrv <- lrv[upper.tri(lrv)]
  
  return(lrv)
}