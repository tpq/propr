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
updateF <- function(propd,
                    moderated = FALSE,
                    ivar = "clr") {
  # Check that active theta is theta_d? propd@active
  if (!propd@active == "theta_d") {
    stop("Make theta_d the active theta.")
  }

  if (moderated) {
    packageCheck("limma")

    # A reference is needed for moderation
    propd@counts # Zeros replaced unless alpha provided...
    use <- index_reference(propd@counts, ivar)

    # Establish data with regard to a reference Z
    if (any(propd@counts == 0)) {
      message("Alert: Building reference set with ivar and counts offset by 1.")
      X <- as.matrix(propd@counts + 1)
    } else{
      message("Alert: Building reference set with ivar and counts.")
      X <- as.matrix(propd@counts)
    }

    logX <- log(X)
    z.set <- logX[, use, drop = FALSE]
    z.geo <- rowMeans(z.set)
    if (any(exp(z.geo) == 0))
      stop("Zeros present in reference set.")
    z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
    z <- exp(z.geo)

    # Fit limma-voom to reference-based data
    message("Alert: Calculating weights with regard to reference.")
    z.sr <-
      t(exp(z.lr) * mean(z)) # scale counts by mean of reference
    design <-
      stats::model.matrix(~ . + 0, data = as.data.frame(propd@group))
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

  } else{
    propd@Fivar <- NA # used by updateCutoffs
    N <- length(propd@group) # population-level metric (i.e., N)
    Fstat <-
      (N - 2) * (1 - propd@results$theta) / propd@results$theta
    theta_mod <- as.numeric(NA)
  }

  propd@results$theta_mod <- theta_mod
  propd@results$Fstat <- Fstat

  # Calculate unadjusted p-value (d1 = K - 1; d2 = N - K)
  K <- length(unique(propd@group))
  N <-
    length(propd@group) + propd@dfz # population-level metric (i.e., N)
  propd@results$Pval <-
    stats::pf(Fstat, K - 1, N - K, lower.tail = FALSE)
  propd@results$FDR <-
    stats::p.adjust(propd@results$Pval, method = "BH")

  return(propd)
}
