#' Calculate Theta
#'
#' Calculate differential proportionality measure, theta.
#'  Used by \code{\link{propd}} to build the \code{@@theta}
#'  slot. A numeric \code{alpha} argument will trigger
#'  the Box-Cox transformation.
#'
#' @inheritParams propd
#' @param lrv A numeric vector. A vector of pre-computed
#'  log-ratio variances. Optional parameter.
#' @param only A character string. The name of the theta
#'  type to return if only calculating one theta type.
#'  Used to make \code{updateCutoffs} faster.
#' @return A \code{data.frame} of \code{theta} values.
#'
#' @export
calculateTheta <- function(counts, group, alpha, lrv = NA, only = "all",
                           weighted = FALSE){

  ct <- as.matrix(counts)
  if(missing(alpha)) alpha <- NA
  if(!is.character(group)) group <- as.character(group)
  if(length(unique(group)) != 2) stop("Please use exactly two unique groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few group labels.")
  if(identical(lrv, NA)){ firstpass <- TRUE
  }else{ firstpass <- FALSE }

  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]
  n1 <- sum(group1)
  n2 <- sum(group2)

  # Calculate weights and lrv modifier
  if(weighted){

    packageCheck("limma")
    design <- matrix(0, nrow = nrow(ct), ncol = 2)
    design[group1, 1] <- 1
    design[group2, 2] <- 1
    v <- limma::voom(t(counts), design = design)
    W <- t(v$weights)
    p1 <- lrvMod(ct[group1,], W[group1,])
    p2 <- lrvMod(ct[group2,], W[group2,])
    p <- lrvMod(ct, W)

  }else{

    W <- ct # not used by lrv() or lrm()
    p1 <- n1 - 1
    p2 <- n2 - 1
    p <- n1 + n2 - 1
  }

  if(is.na(alpha)){

    # Calculate LRV and LRM with provided counts
    if(firstpass) lrv <- lrv(ct, W, weighted)
    lrv1 <- lrv(ct[group1,], W[group1,], weighted)
    lrv2 <- lrv(ct[group2,], W[group2,], weighted)
    lrm1 <- lrm(ct[group1,], W[group1,], weighted)
    lrm2 <- lrm(ct[group2,], W[group2,], weighted)

  }else{

    # Calculate LRV and LRM with 0s present
    if(weighted) stop("No method available for weighted aVLR.")
    if(firstpass){
      if(any(ct == 0)) message("Alert: Approximating LRV in setting of 0s.")
      lrv <- boxRcpp(ct[], alpha)
    }
    lrv1 <- boxRcpp(ct[group1,], alpha)
    lrv2 <- boxRcpp(ct[group2,], alpha)

    # Replace 0s to calculate LRM
    if(firstpass){
      if(any(ct == 0)){
        message("Alert: Replacing 0s to calculate LRM.")
        ct[ct == 0] <- 1 # ct not used again in scope
      }
    }
    lrm1 <- lrm(ct[group1,], W[group1,], weighted)
    lrm2 <- lrm(ct[group2,], W[group2,], weighted)
  }

  # Replace NaN thetas (from VLR = 0) with 1
  lrv0 <- lrv == 0
  replaceNaNs <- any(lrv0)
  if(replaceNaNs){
    if(firstpass) message("Alert: Replacing NaN theta values with 1.")
  }

  # Build all theta types unless only != "all"
  if(only == "all" | only == "theta_d"){

    theta <- (p1 * lrv1 + p2 * lrv2) / (p * lrv)
    if(replaceNaNs) theta[lrv0] <- 1
    if(only == "theta_d") return(theta)
  }

  if(only == "all" | only == "theta_e"){

    theta_e <- 1 - pmax(p1 * lrv1, p2 * lrv2) / (p * lrv)
    if(replaceNaNs) theta_e[lrv0] <- 1
    if(only == "theta_e") return(theta_e)
  }

  if(only == "all" | only == "theta_f"){

    theta_f <- pmax(p1 * lrv1, p2 * lrv2) / (p * lrv)
    if(replaceNaNs) theta_f[lrv0] <- 1
    if(only == "theta_f") return(theta_f)
  }

  F_d <- (n1 + n2 - 2) * (1 - theta) / theta
  labels <- labRcpp(ncol(counts))

  return(
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      "theta" = theta,
      "theta_e" = theta_e,
      "theta_f" = theta_f,
      "lrv" = lrv,
      "lrv1" = lrv1,
      "lrv2" = lrv2,
      "lrm1" = lrm1,
      "lrm2" = lrm2,
      "F_d" = F_d,
      "p1" = p1,
      "p2" = p2,
      "p" = p
    ))
}

#' @rdname propd
#' @export
updateCutoffs <- function(propd, cutoff = seq(.05, .95, .3)){

  # Let NA cutoff skip function
  if(identical(cutoff, NA)) return(propd)

  # Set up FDR cutoff table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
  FDR$cutoff <- cutoff
  p <- ncol(propd@permutes)
  lrv <- propd@theta$lrv

  # Use calculateTheta to permute active theta
  for(k in 1:p){

    numTicks <- progress(k, p, numTicks)

    # Tally k-th thetas that fall below each cutoff
    shuffle <- propd@permutes[, k]
    pkt <- calculateTheta(propd@counts[shuffle, ], propd@group, propd@alpha, lrv,
                          only = propd@active, weighted = propd@weighted)
    for(cut in 1:nrow(FDR)){
      FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + sum(pkt < FDR[cut, "cutoff"])
    }
  }

  # Calculate FDR based on real and permuted tallys
  FDR$randcounts <- FDR$randcounts / p
  for(cut in 1:nrow(FDR)){
    FDR[cut, "truecounts"] <- sum(propd@theta$theta < FDR[cut, "cutoff"])
    FDR[cut, "FDR"] <- FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
  }

  propd@fdr <- FDR
  return(propd)
}
