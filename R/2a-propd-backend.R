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
    NVTX_PUSH("calculate_theta", 0)
    # count matrix
    ct <- as.matrix(counts)

    # define recompute lrv or not
    if (identical(lrv, NA)) {
      firstpass <- TRUE
    } else {
      firstpass <- FALSE
    }

    # Get groups
    NVTX_PUSH("Get groups", 0)
    groups <- lapply(unique(group), function(g) g == group)
    ngrp <- length(unique(group))
    NVTX_POP()

    ##############################################################################
    ### Use weights
    ##############################################################################

    # calculate weights using limma sample weights
    if (weighted) {
      message("Alert: Calculating Limma's reliability weights for samples.")
      packageCheck("limma")

      # use clr-transform of the counts for quality weights from limma:
      NVTX_PUSH("stats::model.matrix", 0)
      design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
      NVTX_POP()

      # limma requires a pseudocount of +1 if zeros are present
      if (any(counts == 0)) {
        X <- as.matrix(counts + 1)
      } else {
        X <- as.matrix(counts)
      }

      # calculate geometric mean
      NVTX_PUSH("geometric mean", 0)
      logX <- log(X)
      z.geo <- rowMeans(logX)
      NVTX_POP()

      # calculate clr-transformed data
      z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
      # scale counts by mean of geometric mean, this rescales data to similar magnitude as original counts
      lz.sr <- t(z.lr + mean(z.geo)) # corresponds to log(z.sr) in updateF function

      # use quality weights from limma:
      NVTX_PUSH("limma::arrayWeights", 0)
      aw <- limma::arrayWeights(lz.sr, design)
      weights <- t(sweep(matrix(1, nrow(lz.sr), ncol(lz.sr)), 2, aw, `*`)) # get the correct dimensions
      NVTX_POP()

      if (nrow(weights) != nrow(counts) | ncol(weights) != ncol(counts)) {
        NVTX_POP() # pop "calculate_theta" before error
        stop("The matrix dimensions of 'weights' must match the matrix dimensions 'counts'.")
      }
    }

    # use weights for lrv modifiers, if provided
    if (weighted) {
      NVTX_PUSH("weighted lrv modifiers", 0)
      W <- weights
      ps <- lapply(groups, function(g) omega(W[g, ]))
      names(ps) <- paste0("p", 1:ngrp)
      p <- omega(W)
      NVTX_POP()
    } else {
      NVTX_PUSH("unweighted lrv modifiers", 0)
      W <- ct
      ps <- lapply(groups, function(g) sum(g) - 1)
      names(ps) <- paste0("p", 1:ngrp)
      p <- length(group) - 1
      NVTX_POP()
    }

    ##############################################################################
    ### Handle zeros
    ##############################################################################

    # Replace zeros if any
    # Logratio and theta cannot be computed with zeros
    if (any(ct == 0) && !is.na(alpha)) {
      NVTX_PUSH("simple_zero_replacement", 0)
      ct <- simple_zero_replacement(ct)
      NVTX_POP()
    }

    ##############################################################################
    ### Calculate logratio variance and differential proportionality theta
    ##############################################################################

    # calculate logratio variance based on shrunk covariance matrix
    if (shrink) {
      NVTX_PUSH("lrv_with_shrinkage", 0)
      if (weighted) {
        NVTX_POP()        # pop "lrv_with_shrinkage"
        NVTX_POP()        # pop "calculate_theta"
        stop("Shrinkage is not available for weighted computation yet.")
      }
      if (firstpass) {
        lrv <- lrv_with_shrinkage(ct)
      }
      lrvs <- lapply(groups, function(g) lrv_with_shrinkage(ct[g, ]))
      names(lrvs) <- paste0("lrv", 1:ngrp)
      NVTX_POP()
    } else {
      # Calculate weighted and/or alpha-transformed LRVs -- W not used if weighted = FALSE
      NVTX_PUSH("lrv", 0)
      if (firstpass) {
        lrv <- lrv(ct, W, weighted, alpha, ct, W)
      }
      lrvs <-
        lapply(groups, function(g)
          lrv(ct[g, ], W[g, ], weighted, alpha, ct, W))
      names(lrvs) <- paste0("lrv", 1:ngrp)
      NVTX_POP()
    }

    # Calculate LRM (using alpha-based LRM if appropriate)
    if (only == "all") {
      NVTX_PUSH("lrm", 0)
      if (firstpass) {
        lrm <- lrm(ct, W, weighted, alpha, ct, W)
      }
      lrms <-
        lapply(groups, function(g)
          lrm(ct[g, ], W[g, ], weighted, alpha, ct, W))
      names(lrms) <- paste0("lrm", 1:ngrp)
      NVTX_POP()
    }

    # Replace NaN thetas (from VLR = 0 or VLR = NaN) with 1
    NVTX_PUSH("replaceNaNs_setup", 0)
    lrv0 <-
      Reduce("|", lapply(lrvs, is.na)) |
      is.na(lrv) | (lrv == 0) # aVLR triggers NaN
    replaceNaNs <- any(lrv0)
    if (replaceNaNs) {
      if (firstpass)
        message("Alert: Replacing NaN theta values with 1.")
    }
    NVTX_POP()

    # Calculate within-group sums-of-squares (used to calculate theta)
    NVTX_PUSH("within_group_SS", 0)
    SS <- lapply(1:ngrp, function(i)
      ps[[i]] * lrvs[[i]])
    NVTX_POP()

    ##############################################################################
    ### Build data.frame of results with computed theta
    ##############################################################################

    # Build all theta types unless only != "all"
    if (only == "all" | only == "theta_d") {
      NVTX_PUSH("theta_d", 0)
      theta <- Reduce(`+`, SS) / (p * lrv)
      if (replaceNaNs)
        theta[lrv0] <- 1
      if (only == "theta_d") {
        NVTX_POP()  # pop "theta_d"
        NVTX_POP()  # pop "calculate_theta"
        return(theta)
      }
      NVTX_POP()    # pop "theta_d"
    }

    if (only == "all" | only == "theta_e") {
      NVTX_PUSH("theta_e", 0)
      theta_e <- 1 - Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_e[lrv0] <- 1
      if (only == "theta_e") {
        NVTX_POP()  # pop "theta_e"
        NVTX_POP()  # pop "calculate_theta"
        return(theta_e)
      }
      NVTX_POP()    # pop "theta_e"
    }

    if (only == "all" | only == "theta_f") {
      NVTX_PUSH("theta_f", 0)
      theta_f <- Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_f[lrv0] <- 1
      if (only == "theta_f") {
        NVTX_POP()  # pop "theta_f"
        NVTX_POP()  # pop "calculate_theta"
        return(theta_f)
      }
      NVTX_POP()    # pop "theta_f"
    }

    if (only == "all" | only == "theta_g") {
      NVTX_PUSH("theta_g", 0)
      theta_g <- Reduce("pmin", SS) / (p * lrv)
      if (replaceNaNs)
        theta_g[lrv0] <- 1
      if (only == "theta_g") {
        NVTX_POP()  # pop "theta_g"
        NVTX_POP()  # pop "calculate_theta"
        return(theta_g)
      }
      NVTX_POP()    # pop "theta_g"
    }

    NVTX_PUSH("build_result", 0)
    labels <- labRcpp(ncol(counts))
    res <- data.frame(
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
    NVTX_POP()      # pop "build_result"
    NVTX_POP()      # pop "calculate_theta"
    return(res)
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
  NVTX_PUSH("lrv_with_shrinkage", 0)
  
  # compute covariance matrix on the log data
  NVTX_PUSH("compute_covariance", 0)
  if (shrink) {
    NVTX_PUSH("corpcor::cov.shrink", 0)
    cov_matrix <- corpcor::cov.shrink(log(ct))
    NVTX_POP()
  } else {
    NVTX_PUSH("stats::cov", 0)
    cov_matrix <- stats::cov(log(ct))
    NVTX_POP()
  }
  NVTX_POP()  # compute_covariance

  # convert shrunk covariance matrix to shrunk logratio variance matrix
  NVTX_PUSH("build_lrv_matrix", 0)
  diag_matrix <- diag(cov_matrix)
  outer_sum <- outer(diag_matrix, diag_matrix, "+")
  lrv <- outer_sum - 2 * cov_matrix
  NVTX_POP()

  # get upper triangle so it is coherent with propr:::lrv
  NVTX_PUSH("upper_triangle", 0)
  lrv <- lrv[upper.tri(lrv)]
  NVTX_POP()
  
  NVTX_POP()  # lrv_with_shrinkage
  return(lrv)
}
