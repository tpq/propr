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
#' @param weights A weight matrix.
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
           weights = as.matrix(NA),
           shrink = FALSE) {

    NVTX_PUSH("calculate_theta_inside", 1)
    
    # count matrix
    NVTX_PUSH("count matrix", 2)
    ct <- as.matrix(counts)
    NVTX_POP()

    # define recompute lrv or not
    NVTX_PUSH("parameter_setup", 2)
    if (identical(lrv, NA)) {
      firstpass <- TRUE
    } else{
      firstpass <- FALSE
    }
    NVTX_POP()

    # Get groups
    NVTX_PUSH("get_groups", 2)
    groups <- lapply(unique(group), function(g) g == group)
    ngrp   <- length(unique(group))
    NVTX_POP()

    ##############################################################################
    ### Use weights
    ##############################################################################

    if (weighted) {
      # calculate weights using limma
      if (is.na(weights[1,1])) {
        message("Alert: Calculating limma-based weights.")
        packageCheck("limma")
        NVTX_PUSH("build_design_and_voom", 2)
        design  <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
        v       <- limma::voom(t(counts), design = design)
        weights <- t(v$weights)
        NVTX_POP()
      }
      
      NVTX_PUSH("weight_matrix_validation", 2)
      if (nrow(weights) != nrow(counts) | ncol(weights) != ncol(counts)) {
        stop("The matrix dimensions of 'weights' must match the matrix dimensions 'counts'.")
      }
      NVTX_POP()
    }

    # use weights for lrv modifiers, if provided
    if (weighted) {
      W <- weights
      message("omega group (lrv modifiers )")
      #batch apply omega to the weights
      NVTX_PUSH("omega group (lrv modifiers )", 2)
      ps <- lapply(groups, function(g) omega(W[g,]))
      NVTX_POP()
      message("DONE: omega group (lrv modifiers )")

      names(ps) <- paste0("p", 1:ngrp)
      message("omega single (lrv modifiers )")
      NVTX_PUSH("omega single (lrv modifiers )", 2)
      p <- omega(W)
      NVTX_POP()
      message("DONE: omega single (lrv modifiers )")

    } else {
      NVTX_PUSH("compute_W_and_ps", 2)
      W <- ct
      ps <- lapply(groups, function(g) sum(g) - 1)
      names(ps) <- paste0("p", 1:ngrp)
      p <- length(group) - 1
      NVTX_POP()
    }

    ##############################################################################
    ### Calculate logratio variance and differential proportionality theta
    ##############################################################################

    # calculate logratio variance based on shrunk covariance matrix
    if (shrink) {
      NVTX_PUSH("shrinkage_validation", 2)
      if (weighted) {
        stop("Shrinkage is not available for weighted computation yet.")
      }
      NVTX_POP()
      
      if (firstpass) {
          NVTX_PUSH("lrv_with_shrinkage_single (firstpass)", 2)
            lrv <- lrv_with_shrinkage(ct)
          NVTX_POP()
      }

      NVTX_PUSH("lrv_with_shrinkage_group", 2)
        lrvs <-  lapply(groups, function(g) lrv_with_shrinkage(ct[g,]))
      NVTX_POP()

      NVTX_PUSH("name_lrvs_shrinkage", 2)
      names(lrvs) <- paste0("lrv", 1:ngrp)
      NVTX_POP()
    } else {
        # Calculate weighted and/or alpha-transformed LRVs -- W not used if weighted = FALSE
        if (firstpass) {
          NVTX_PUSH("lrv_single (first_pass)", 2)
          lrv <- lrv(ct, W, weighted, alpha, ct, W)
          NVTX_POP()
        }
        NVTX_PUSH("lrv_group", 2)
        lrvs <- lapply(groups, function(g) {
          lrv(ct[g,], W[g,], weighted, alpha, ct, W)
        })
        NVTX_POP()
        
        NVTX_PUSH("name_lrvs", 2)
        names(lrvs) <- paste0("lrv", 1:ngrp)
        NVTX_POP()
    }

    # Calculate LRM (using alpha-based LRM if appropriate)
    if (only == "all") {
      if (firstpass) {
        NVTX_PUSH("lrm_single (first_pass)", 2)
        lrm <- lrm(ct, W, weighted, alpha, ct, W)
        NVTX_POP()
      }
      NVTX_PUSH("lrm_group", 2)
        lrms <- lapply(groups, function(g) {
          lrm(ct[g,], W[g,], weighted, alpha, ct, W)
        })
      NVTX_POP()
      
      NVTX_PUSH("name_lrms", 2)
      names(lrms) <- paste0("lrm", 1:ngrp)
      NVTX_POP()
    }

    # Replace NaN thetas (from VLR = 0 or VLR = NaN) with 1
    NVTX_PUSH("nan_replacement_setup", 2)
    lrv0 <- Reduce("|", lapply(lrvs, is.na)) | is.na(lrv) | (lrv == 0) # aVLR triggers NaN
    replaceNaNs <- any(lrv0)
    if (replaceNaNs) {
      if (firstpass) message("Alert: Replacing NaN theta values with 1.")
    }
    NVTX_POP()

    # Calculate within-group sums-of-squares (used to calculate theta)
    NVTX_PUSH("sum_of_squares_calculation", 2)
    SS <- lapply(1:ngrp, function(i)
      ps[[i]] * lrvs[[i]])
    NVTX_POP()

    ##############################################################################
    ### Build data.frame of results with computed theta
    ##############################################################################

    # Build all theta types unless only != "all"
    if (only == "all" | only == "theta_d") {
      NVTX_PUSH("calculate_theta_d", 2)
      theta <- Reduce(`+`, SS) / (p * lrv)
      if (replaceNaNs)
        theta[lrv0] <- 1
      NVTX_POP()
      
      if (only == "theta_d") {
        NVTX_POP()  # Close main function scope
        return(theta)
      }
    }

    if (only == "all" | only == "theta_e") {
      NVTX_PUSH("calculate_theta_e", 2)
      theta_e <- 1 - Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_e[lrv0] <- 1
      NVTX_POP()
      
      if (only == "theta_e") {
        NVTX_POP()  # Close main function scope
        return(theta_e)
      }
    }

    if (only == "all" | only == "theta_f") {
      NVTX_PUSH("calculate_theta_f", 2)
      theta_f <- Reduce("pmax", SS) / (p * lrv)
      if (replaceNaNs)
        theta_f[lrv0] <- 1
      NVTX_POP()
      
      if (only == "theta_f") {
        NVTX_POP()  # Close main function scope
        return(theta_f)
      }
    }

    if (only == "all" | only == "theta_g") {
      NVTX_PUSH("calculate_theta_g", 2)
      theta_g <- Reduce("pmin", SS) / (p * lrv)
      if (replaceNaNs)
        theta_g[lrv0] <- 1
      NVTX_POP()
    }

    NVTX_PUSH("build_final_dataframe", 2)
    labels <- labRcpp(ncol(counts))
    result <- data.frame(
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
    NVTX_POP()

    NVTX_POP()
    return(result)
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
  
  NVTX_PUSH("lrv_with_shrinkage_function", 1)
  
  # compute covariance matrix on the log data
  NVTX_PUSH("compute_covariance_matrix", 2)
  if (shrink) {
    cov_matrix <- corpcor::cov.shrink(log(ct))
  } else {
    cov_matrix <- stats::cov(log(ct))
  }
  NVTX_POP()

  # convert shrunked covariance matrix to shrunked logratio variance matrix
  NVTX_PUSH("convert_to_lrv_matrix", 2)
  diag_matrix <- diag(cov_matrix)
  outer_sum <- outer(diag_matrix, diag_matrix, "+")
  lrv <- outer_sum - 2 * cov_matrix
  NVTX_POP()

  # it is symmetric
  # we get the upper triangle, so that it is coherent with the function propr:::lrv
  NVTX_PUSH("extract_upper_triangle", 2)
  lrv <- lrv[upper.tri(lrv)]
  NVTX_POP()
  
  NVTX_POP()
  return(lrv)
}