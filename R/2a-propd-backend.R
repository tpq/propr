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
#'  should be performed. If `TRUE`, the function will use limma-based weights
#'  for the calculations.
#' @param weights A matrix of limma-based weights. This parameter is optional
#'  and used only if `weighted = TRUE`.
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
           weights = as.matrix(NA)) {
    ct <- as.matrix(counts)
    if (identical(lrv, NA)) {
      firstpass <- TRUE
    } else{
      firstpass <- FALSE
    }

    # Get groups
    groups <- lapply(unique(group), function(g)
      g == group)
    ngrp <- length(unique(group))

    # Calculate weights and lrv modifier
    if (weighted) {
      # Do not delete -- this is used by updateCutoffs.propd()
      if (is.na(weights[1, 1])) {
        message("Alert: Calculating limma-based weights.")
        packageCheck("limma")
        design <-
          stats::model.matrix(~ . + 0, data = as.data.frame(group))
        v <- limma::voom(t(counts), design = design)
        weights <- t(v$weights)
      }

      W <- weights
      ps <- lapply(groups, function(g)
        omega(W[g,]))
      names(ps) <- paste0("p", 1:ngrp)
      p <- omega(W)

    } else{
      W <- ct
      ps <- lapply(groups, function(g)
        sum(g) - 1)
      names(ps) <- paste0("p", 1:ngrp)
      p <- length(group) - 1
    }

    # Calculate weighted and/or alpha-transformed LRVs -- W not used if weighted = FALSE
    if (firstpass)
      lrv <- lrv(ct, W, weighted, alpha, ct, W)
    lrvs <-
      lapply(groups, function(g)
        lrv(ct[g,], W[g,], weighted, alpha, ct, W))
    names(lrvs) <- paste0("lrv", 1:ngrp)

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

#' Parse propd theta values to a matrix
#'
#' This function returns a symmetric matrix of \code{propd} values.
#' Note that the active theta values are used for the matrix. The 
#' diagonal is set to NA.
#'
#' @inheritParams getResults
#'
#' @return A symmetric matrix.
#'
#' @export
get_theta_matrix <- function(propd){

  if (!inherits(propd, "propd"))
      stop("Please provide a 'propd' object.")
  if(propd@results$Partner[1] != 2 | propd@results$Pair[1] != 1){
    stop("Unexpected sorting of results slot.")
  }

  # Convert to matrix
  mat <- half2mat(propd@results[,'theta'])
  rownames(mat) <- colnames(propd@counts)
  colnames(mat) <- colnames(propd@counts)

  # TODO check if this is the correct way to handle this
  # and why before is 0
  diag(mat) <- NA

  return(mat)
}