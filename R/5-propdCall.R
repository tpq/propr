#' Calculate Theta
#'
#' Calculate differential proportionality measure, theta.
#'  Used by \code{\link{propd}} to build the \code{@@results}
#'  slot. A numeric \code{alpha} argument will trigger
#'  the Box-Cox transformation.
#'
#' @inheritParams all
#' @param lrv A numeric vector. A vector of pre-computed
#'  log-ratio variances. Optional parameter.
#' @param only A character string. The name of the theta
#'  type to return if only calculating one theta type.
#'  Used to make \code{updateCutoffs} faster.
#' @param weights A matrix. Pre-computed \code{limma}-based
#'  weights. Optional parameter.
#'
#' @return A data.frame of theta values if \code{only = "all"}.
#'  Otherwise, this function returns a numeric vector.
calculateTheta <- function(counts, group, alpha = NA, lrv = NA, only = "all",
                           weighted = FALSE, weights = as.matrix(NA)){

  ct <- as.matrix(counts)
  if(identical(lrv, NA)){ firstpass <- TRUE
  }else{ firstpass <- FALSE }

  # Get groups
  groups <- lapply(unique(group), function(g) g == group)
  ngrp <- length(unique(group))

  # Calculate weights and lrv modifier
  if(weighted){

    # Do not delete -- this is used by updateCutoffs.propd()
    if(is.na(weights[1,1])){
      message("Alert: Calculating limma-based weights.")
      packageCheck("limma")
      design <- stats::model.matrix(~.+0, data = as.data.frame(group))
      v <- limma::voom(t(counts), design = design)
      weights <- t(v$weights)
    }

    W <- weights
    ps <- lapply(groups, function(g) omega(W[g,]))
    names(ps) <- paste0("p", 1:ngrp)
    p <- omega(W)

  }else{

    W <- ct
    ps <- lapply(groups, function(g) sum(g)-1)
    names(ps) <- paste0("p", 1:ngrp)
    p <- length(group) - 1
  }

  # Calculate weighted and/or alpha-transformed LRVs -- W not used if weighted = FALSE
  if(firstpass) lrv <- lrv(ct, W, weighted, alpha, ct, W)
  lrvs <- lapply(groups, function(g) lrv(ct[g,], W[g,], weighted, alpha, ct, W))
  names(lrvs) <- paste0("lrv", 1:ngrp)

  # Calculate LRM (using alpha-based LRM if appropriate)
  if(only == "all"){
    if(firstpass) lrm <- lrm(ct, W, weighted, alpha, ct, W)
    lrms <- lapply(groups, function(g) lrm(ct[g,], W[g,], weighted, alpha, ct, W))
    names(lrms) <- paste0("lrm", 1:ngrp)
  }

  # Replace NaN thetas (from VLR = 0 or VLR = NaN) with 1
  lrv0 <- Reduce("|", lapply(lrvs, is.na)) | is.na(lrv) | (lrv == 0) # aVLR triggers NaN
  replaceNaNs <- any(lrv0)
  if(replaceNaNs){
    if(firstpass) message("Alert: Replacing NaN theta values with 1.")
  }

  # Calculate within-group sums-of-squares (used to calculate theta)
  SS <- lapply(1:ngrp, function(i) ps[[i]] * lrvs[[i]])

  # Build all theta types unless only != "all"
  if(only == "all" | only == "theta_d"){

    theta <- Reduce(`+`, SS) / (p * lrv)
    if(replaceNaNs) theta[lrv0] <- 1
    if(only == "theta_d") return(theta)
  }

  if(only == "all" | only == "theta_e"){

    theta_e <- 1 - Reduce("pmax", SS) / (p * lrv)
    if(replaceNaNs) theta_e[lrv0] <- 1
    if(only == "theta_e") return(theta_e)
  }

  if(only == "all" | only == "theta_f"){

    theta_f <- Reduce("pmax", SS) / (p * lrv)
    if(replaceNaNs) theta_f[lrv0] <- 1
    if(only == "theta_f") return(theta_f)
  }

  if(only == "all" | only == "theta_g"){

    theta_g <- Reduce("pmin", SS) / (p * lrv)
    if(replaceNaNs) theta_g[lrv0] <- 1
    if(only == "theta_g") return(theta_g)
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
    ))
}

#' @rdname propd
#' @section Functions:
#' \code{updateCutoffs:}
#'  Use the \code{propd} object to permute theta across a
#'  number of theta cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#' @export
updateCutoffs.propd <- function(object, cutoff = seq(.05, .95, .3)){

  if(identical(object@permutes, data.frame())) stop("Permutation testing is disabled.")

  # Let NA cutoff skip function
  if(identical(cutoff, NA)) return(object)

  # Set up FDR cutoff table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
  FDR$cutoff <- cutoff
  p <- ncol(object@permutes)
  lrv <- object@results$lrv

  # Use calculateTheta to permute active theta
  for(k in 1:p){

    numTicks <- progress(k, p, numTicks)

    # Tally k-th thetas that fall below each cutoff
    shuffle <- object@permutes[, k]

    if(object@active == "theta_mod"){

      # Calculate theta_mod with updateF (using i-th permuted object)
      if(is.na(object@Fivar)) stop("Please re-run 'updateF' with 'moderation = TRUE'.")
      propdi <- suppressMessages(
        propd(object@counts[shuffle, ], group = object@group, alpha = object@alpha, p = 0,
              weighted = object@weighted))
      propdi <- suppressMessages(
        updateF(propdi, moderated = TRUE, ivar = object@Fivar))
      pkt <- propdi@results$theta_mod

    }else{

      # Calculate all other thetas directly (using calculateTheta)
      pkt <- suppressMessages(
        calculateTheta(object@counts[shuffle, ], object@group, object@alpha, lrv,
                       only = object@active, weighted = object@weighted))
    }

    # Find number of permuted theta less than cutoff
    for(cut in 1:nrow(FDR)){ # randcounts as cumsum
      FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + sum(pkt < FDR[cut, "cutoff"])
    }
  }

  # Calculate FDR based on real and permuted tallys
  FDR$randcounts <- FDR$randcounts / p # randcounts as mean
  for(cut in 1:nrow(FDR)){
    FDR[cut, "truecounts"] <- sum(object@results$theta < FDR[cut, "cutoff"])
    FDR[cut, "FDR"] <- FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
  }

  # Initialize @fdr
  object@fdr <- FDR

  return(object)
}

#' @rdname propd
#' @section Functions:
#' \code{updateF:}
#'  Use the \code{propd} object to calculate the F-statistic
#'  from theta as described in the Erb et al. 2017 manuscript
#'  on differential proportionality. Optionally calculates a
#'  moderated F-statistic using the limma-voom method. Supports
#'  weighted and alpha transformed theta values.
#' @export
updateF <- function(propd, moderated = FALSE, ivar = "clr"){

  # Check that active theta is theta_d? propd@active
  if(!propd@active == "theta_d"){
    stop("Make theta_d the active theta.")
  }

  if(moderated){

    packageCheck("limma")

    # A reference is needed for moderation
    propd@counts # Zeros replaced unless alpha provided...
    use <- ivar2index(propd@counts, ivar)

    # Establish data with regard to a reference Z
    if(any(propd@counts == 0)){
      message("Alert: Building reference set with ivar and counts offset by 1.")
      X <- as.matrix(propd@counts + 1)
    }else{
      message("Alert: Building reference set with ivar and counts.")
      X <- as.matrix(propd@counts)
    }

    logX <- log(X)
    z.set <- logX[, use, drop = FALSE]
    z.geo <- rowMeans(z.set)
    if(any(exp(z.geo) == 0)) stop("Zeros present in reference set.")
    z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
    z <- exp(z.geo)

    # Fit limma-voom to reference-based data
    message("Alert: Calculating weights with regard to reference.")
    z.sr <- t(exp(z.lr) * mean(z)) # scale counts by mean of reference
    design <- stats::model.matrix(~.+0, data = as.data.frame(propd@group))
    v <- limma::voom(z.sr, design = design)
    param <- limma::lmFit(v, design)
    param <- limma::eBayes(param)
    z.df <- param$df.prior
    propd@dfz <- param$df.prior
    z.s2 <- param$s2.prior

    # Calculate simple moderation term based only on LRV
    mod <- z.df * z.s2 / propd@results$lrv

    # Moderate F-statistic
    propd@Fivar <- ivar # used by updateCutoffs
    N <- length(propd@group) # population-level metric (i.e., N)
    Fprime <- (1 - propd@results$theta) * (N + z.df) /
      ((N) * propd@results$theta + mod)
    Fstat <- (N + z.df - 2) * Fprime
    theta_mod <- 1 / (1 + Fprime)

  }else{

    propd@Fivar <- NA # used by updateCutoffs
    N <- length(propd@group) # population-level metric (i.e., N)
    Fstat <- (N - 2) * (1 - propd@results$theta) / propd@results$theta
    theta_mod <- as.numeric(NA)
  }

  propd@results$theta_mod <- theta_mod
  propd@results$Fstat <- Fstat

  # Calculate unadjusted p-value (d1 = K - 1; d2 = N - K)
  K <- length(unique(propd@group))
  N <- length(propd@group) + propd@dfz # population-level metric (i.e., N)
  propd@results$Pval <- pf(Fstat, K - 1, N - K, lower.tail = FALSE)
  propd@results$FDR <- p.adjust(propd@results$Pval, method = "BH")

  return(propd)
}

#' Calculate a theta Cutoff
#'
#' This function uses the F distribution to calculate a cutoff of
#'  theta for a p-value given by the \code{pval} argument.
#'
#' If the argument \code{fdr = TRUE}, this function returns the
#'  empiric cutoff that corresponds to the FDR-adjusted p-value
#'  stored in the \code{@@results$FDR} slot.
#'
#' @inheritParams all
#' @param pval A p-value at which to calculate a theta cutoff.
#' @param fdr A boolean. Toggles whether to calculate the theta
#'  cutoff for an FDR-adjusted p-value.
#'
#' @return A cutoff of theta from [0, 1].
#'
#' @export
qtheta <- function(propd, pval = 0.05, fdr = FALSE){

  if(!"Fstat" %in% colnames(propd@results)){
    stop("Please run updateF() on propd object before calling qtheta.")
  }

  if(pval < 0 | pval > 1){
    stop("Provide a p-value cutoff from [0, 1].")
  }

  if(fdr){

    message("Alert: Returning an empiric cutoff based on the $FDR slot.")
    index <- propd@results$FDR < pval
    if(any(index)){
      cutoff <- max(propd@results$theta[index])
    }else{
      stop("No pairs below p-value.")
    }

  }else{

    # Compute based on theory
    K <- length(unique(propd@group))
    N <- length(propd@group) + propd@dfz # population-level metric (i.e., N)
    Q <- qf(pval, K - 1, N - K, lower.tail = FALSE)
    # # Fstat <- (N - 2) * (1 - propd@theta$theta) / propd@theta$theta
    # # Q = Fstat
    # # Q = (N-2) * (1-theta) / theta
    # # Q / (N-2) = (1/theta) - 1
    # # 1/theta = Q / (N-2) + 1 = Q(N-2)/(N-2)
    # # theta = (N-2)/(Q+(N-2))
    cutoff <- (N-2)/(Q+(N-2))
  }

  return(cutoff)
}
