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
#'  type to calculate if only calculating one theta type.
#'  Used to make \code{updateCutoffs} faster.
#' @return A \code{data.frame} of \code{theta} values.
#'
#' @export
calculateTheta <- function(counts, group, alpha, lrv, only = "all"){

  ct <- as.matrix(counts)
  if(missing(alpha)) alpha <- NA
  if(!is.character(group)) group <- as.character(group)
  if(length(unique(group)) != 2) stop("Please use exactly two unique groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few group labels.")
  if(missing(lrv)){ firstpass <- TRUE
  }else{ firstpass <- FALSE }

  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]
  n1 <- sum(group1)
  n2 <- sum(group2)

  if(is.na(alpha)){

    if(firstpass) lrv <- lltRcpp(vlrRcpp(ct[]))
    lrv1 <- lltRcpp(vlrRcpp(ct[group1,]))
    lrv2 <- lltRcpp(vlrRcpp(ct[group2,]))

  }else{

    if(firstpass) lrv <- boxRcpp(ct[], alpha)
    lrv1 <- boxRcpp(ct[group1,], alpha)
    lrv2 <- boxRcpp(ct[group2,], alpha)
  }

  # Replace NaN thetas (from VLR = 0) with 1
  lrv0 <- lrv == 0
  replaceNaNs <- any(lrv0)
  if(replaceNaNs){
    if(firstpass) message("Alert: Replacing NaN theta values with 1.")
  }

  # Build all theta types unless only != "all"
  if(only == "all" | only == "d"){

    theta <- ((n1-1) * lrv1 + (n2-1) * lrv2) / ((n1+n2-1) * lrv)
    if(replaceNaNs) theta[lrv0] <- 1
    if(only == "d") return(theta)
  }

  if(only == "all" | only == "e"){

    theta_e <- 1 - pmax((n1-1) * lrv1, (n2-1) * lrv2) / ((n1+n2-1) * lrv)
    if(replaceNaNs) theta_e[lrv0] <- 1
    if(only == "e") return(theta_e)
  }

  if(only == "all" | only == "f"){

    theta_f <- 1 - theta_e
    if(replaceNaNs) theta_f[lrv0] <- 1
    if(only == "f") return(theta_f)
  }

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
      "lrv2" = lrv2
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

  # Tally permuted thetas that fall below each cutoff
  i <- which("theta" == colnames(propd@theta))
  if(i == 3){
    cat("Permuting disjointed proportionality (theta_d):\n")
    onlyTheta <- "d"
  }else if(i == 4){
    if(i == 4) cat("Permuting emergent proportionality (theta_e):\n")
    onlyTheta <- "e"
  }else if(i == 5){
    if(i == 5) cat("Permuting fettered proportionality (theta_f):\n")
    onlyTheta <- "f"
  }else{
    stop("No 'updateCutoffs' method in place for active theta.")
  }

  # Use calculateTheta to permute theta
  for(k in 1:p){

    numTicks <- progress(k, p, numTicks)

    # Tally k-th thetas that fall below each cutoff
    shuffle <- propd@permutes[, k]
    pkt <- calculateTheta(propd@counts[shuffle, ], propd@group, propd@alpha, lrv, only = onlyTheta)
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
