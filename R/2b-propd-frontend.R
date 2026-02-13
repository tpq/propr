#' @rdname propd
#' @section Functions:
#' \code{setDisjointed:}
#'  A wrapper for \code{setActive(propd, what = "theta_d")}.
#' @export
setDisjointed <-
  function(propd) {
    setActive(propd, what = "theta_d")
  }

#' @rdname propd
#' @section Functions:
#' \code{setEmergent:}
#'  A wrapper for \code{setActive(propd, what = "theta_e")}.
#' @export
setEmergent <-
  function(propd) {
    if (!all(table(propd@group) == table(propd@group)[1])) {
      warning("Emergent proportionality not yet validated for unequal group sizes.")
    }

    setActive(propd, what = "theta_e")
  }

#' @rdname propd
#' @param what A character string. The theta type to set active.
#' @section Functions:
#' \code{setActive:}
#'  Build analyses and figures using a specific theta type. For
#'  example, set \code{what = "theta_d"} to analyze disjointed
#'  proportionality and \code{what = "theta_e"} to analyze
#'  emergent proportionality.
#' @export
setActive <-
  function(propd, what = "theta_d") {
    if (!inherits(propd, "propd"))
      stop("Please provide a 'propd' object.")
    if (!any(what == colnames(propd@results)) &
        what != propd@active) {
      stop("Provided theta type not recognized.")
    }

    # Rename old active theta type
    i <- which(colnames(propd@results) == "theta")
    colnames(propd@results)[i] <- propd@active

    # Name new active theta type
    i <- which(colnames(propd@results) == what)
    colnames(propd@results)[i] <- "theta"

    propd@active <- what
    message("Alert: Update FDR or F-stat manually.")

    return(propd)
  }

#' @rdname propd
#' @param moderated For \code{updateF}, a boolean. Toggles
#'  whether to calculate a moderated F-statistic.
#' @param ivar See \code{propr} method.
#' @section Functions:
#' \code{updateF:}
#'  Use the \code{propd} object to calculate the F-statistic
#'  from theta as described in the Erb et al. 2017 manuscript
#'  on differential proportionality. Optionally calculates a
#'  moderated F-statistic using the limma-voom method. Supports
#'  weighted and alpha transformed theta values.
#' @export
updateF <- function(propd, moderated = TRUE, ivar = "clr") {
  # Validate inputs
  if (!propd@active == "theta_d") {
    stop("Make theta_d the active theta.")
  }
  
  if (moderated) {
    message("Alert: Calculating moderated F-statistics using voom approach, with trend = ", propd@weighted)

    # use limma package
    packageCheck("limma")
    
    # Set ivar for updateCutoffs
    propd@Fivar <- ivar

    # Prepare reference data
    ref_data <- prepare_reference_data(propd@counts, ivar)

    # Fit limma model to get degrees of freedom
    limma_params <- fit_limma_model(
      ref_data$z.lr, 
      ref_data$z, 
      propd@group, 
      propd@weighted)
      
    z.df <- limma_params$df.prior
    propd@dfz <- z.df
    
    if (propd@weighted) {
      # Calculate variance estimates using new binning approach
      S2 <- calculate_variance_estimates(propd, ref_data$logX)
    } else {
      S2 <- limma_params$s2.prior
    }
      
    # Calculate moderation term
    mod <- z.df * S2 / propd@results$lrv
      
    # Calculate F-statistics
    f_results <- calculate_f_statistics(propd, mod, z.df, moderated = TRUE)

  } else {
    message("Alert: Calculating F-statistics without moderation.")
    propd@Fivar <- NA
    f_results <- calculate_f_statistics(propd, NULL, NULL, moderated = FALSE)
  }
  
  # Store results
  propd@results$theta_mod <- f_results$theta_mod
  propd@results$Fstat <- f_results$Fstat
  
  # Calculate p-values and FDR
  p_results <- calculate_p_values(propd, f_results$Fstat)
  propd@results$Pval <- p_results$Pval
  propd@results$FDR <- p_results$FDR
  
  return(propd)
}


# Helper function to prepare reference data
prepare_reference_data <- function(counts, ivar) {
  use <- index_reference(counts, ivar)
  
  # Handle zeros by adding offset
  if (any(counts == 0)) {
    message("Alert: Building reference set with ivar and counts offset by 1.")
    X <- as.matrix(counts + 1)
  } else {
    message("Alert: Building reference set with ivar and counts.")
    X <- as.matrix(counts)
  }
  
  # Transform to log scale
  logX <- log(X)
  z.set <- logX[, use, drop = FALSE]
  z.geo <- rowMeans(z.set)
  
  if (any(exp(z.geo) == 0)) {
    stop("Zeros present in reference set.")
  }
  
  z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
  z <- exp(z.geo)
  
  return(list(
    logX = logX,
    z.lr = z.lr,
    z = z,
    use = use
  ))
}

# Helper function to calculate variance estimates using binning approach
calculate_variance_estimates <- function(propd, logX) {
  # Pooled log-ratio variance within groups
  plrv <- (propd@results$p1 * propd@results$lrv1 + 
           propd@results$p2 * propd@results$lrv2) / propd@results$p
  
  # Raw count log geometric averages
  avs <- apply(logX, 2, mean)
  
  # Calculate pairwise log geometric mean
  Z <- log(exp(avs) %*% t(exp(avs)))
  add <- Z[upper.tri(Z)]
  
  # Create bins for variance estimation
  nbins <- min(20, ncol(logX))
  As <- sort(add)
  B <- round(length(add) / nbins)
  
  # Define bin boundaries
  se <- rep(0, nbins + 1)
  se[1] <- As[1] - 0.1  # Slightly smaller than min for > comparison
  for (i in 1:(nbins - 1)) {
    se[i + 1] <- As[B * i]
  }
  se[nbins + 1] <- As[length(As)]
  
  # Assign pairs to bins
  L <- list()
  for (i in 1:(length(se) - 1)) {
    L[[i]] <- which(add > se[i] & add <= se[i + 1])
  }
  
  # Calculate smoothed variance estimates
  S <- rep(0, length(L))
  for (i in 1:length(L)) {
    S[i] <- mean(plrv[L[[i]]]^(1/4))  # Similar to limma trend
  }
  
  # Map variance estimates back to original pairs
  S2 <- rep(0, length(add))
  for (i in 1:length(L)) {
    S2[L[[i]]] <- S[i]^4
  }
  
  return(S2)
}

# Helper function to fit limma model and extract parameters
fit_limma_model <- function(z.lr, z, group, moderated_trend = FALSE) {
  if (moderated_trend) {
    message("Alert: Calculating prior degrees of freedom with regard to reference.")
  } else {
    message("Alert: Calculating weights with regard to reference.")
  }
  
  # Scale counts by mean of reference
  z.sr <- t(exp(z.lr) * mean(z))
  
  # Create design matrix
  design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
  
  # Fit limma-voom model
  v <- limma::voom(z.sr, design = design)
  param <- limma::lmFit(v, design)
  param <- limma::eBayes(param, trend = moderated_trend)
  
  return(list(
    df.prior = param$df.prior,
    s2.prior = param$s2.prior
  ))
}

# Helper function to calculate F-statistics
calculate_f_statistics <- function(propd, mod, z.df, moderated = TRUE) {
  N <- length(propd@group)
  
  if (moderated) {
    Fprime <- (1 - propd@results$theta) * (N + z.df) /
              ((N) * propd@results$theta + mod)
    Fstat <- (N + z.df - 2) * Fprime
    theta_mod <- 1 / (1 + Fprime)
  } else {
    Fstat <- (N - 2) * (1 - propd@results$theta) / propd@results$theta
    theta_mod <- as.numeric(NA)
  }
  
  return(list(
    Fstat = Fstat,
    theta_mod = theta_mod
  ))
}

# Helper function to calculate p-values and FDR
calculate_p_values <- function(propd, Fstat) {
  K <- length(unique(propd@group))
  N <- length(propd@group) + propd@dfz
  
  Pval <- stats::pf(Fstat, K - 1, N - K, lower.tail = FALSE)
  FDR <- stats::p.adjust(Pval, method = "BH")
  
  return(list(Pval = Pval, FDR = FDR))
}