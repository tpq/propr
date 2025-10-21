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
#' @param moderated_trend For \code{updateF}, a boolean. Toggles
#'  whether to incorporate mean-variance trend in the moderated 
#'  F-statistic. Default is FALSE.
#' @param ivar See \code{propr} method.
#' @section Functions:
#' \code{updateF:}
#'  Use the \code{propd} object to calculate the F-statistic
#'  from theta as described in the Erb et al. 2017 manuscript
#'  on differential proportionality. Optionally calculates a
#'  moderated F-statistic using the limma-voom method. Supports
#'  weighted and alpha transformed theta values.
#' @export
#check consequences for downward compatibility when replacing moderated=FALSE by moderated=TRUE

updateF <- function(propd,
                    moderated = TRUE,
                    ivar = "clr") {
  # Check that active theta is theta_d? propd@active
  if (!propd@active == "theta_d") {
    stop("Make theta_d the active theta.")
  }

  if (moderated) {
    propr:::packageCheck("limma")

    # A reference is needed for moderation
    propd@counts # Zeros replaced unless alpha provided...
    use <- propr:::index_reference(propd@counts, ivar)

    # Establish data with regard to a reference Z
    if (any(propd@counts == 0)) {
      message("Alert: Building reference set with ivar and counts offset by 1.")
      X <- as.matrix(propd@counts + 1)  # TODO this should be checked, if we still want to do this
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
    message("Alert: Calculating prior degrees of freedom with regard to reference.") #voom weights no longer needed
    z.sr <-
      t(exp(z.lr) * mean(z)) # scale counts by mean of reference
    design <-
      stats::model.matrix(~ . + 0, data = as.data.frame(propd@group))
    v <- limma::voom(z.sr, design = design)
    param <- limma::lmFit(v, design) #check if weights should be used
    param <- limma::eBayes(param,trend=TRUE) #we here fit limma trend genewise with respect to reference, only needed for z.df, not z.s2
    z.df <- param$df.prior
    propd@dfz <- param$df.prior #same as z.df
    #z.s2 <- param$s2.prior #this is no longer needed!

    #new way of calculating z.s2 directly from logratios:
    #pooled lrv within groups:
	
    plrv=(propd@results$p1*propd@results$lrv1+propd@results$p2*propd@results$lrv2)/propd@results$p
    #raw count log geometric averages:
    avs=apply(logX,2,mean)
    #Calculate pairwise log geometric mean adding both genewise avs
    #(do this in the exponent for vectorization of calculations):
    Z=log(exp(avs)%*%t(exp(avs)))
    #write as vector in same format as pairwise results:
    add=Z[upper.tri(Z)]
    nbins=min(20,ncol(logX))
    As=sort(add)
    B=round(length(add)/nbins)
    se=rep(0,nbins+1) #sequence of bin boundaries for entries in As, 20 is arbitrary but works
    se[1]=As[1]-0.1 #smaller than min(As) so we can use > later 
    for (i in 1:(nbins-1)){
	se[i+1]=As[B*i]
    }
    se[nbins+1]=As[length(As)] #max(As)
    L=list()
    for (i in 1:(length(se)-1)){
	L[[i]]=which(add>se[i] & add<=se[i+1])
    }
    S=rep(0,length(L))
    for (i in 1:length(L)){
	S[i]=mean(plrv[L[[i]]]^(1/4)) #similar to limma trend use sqrt(sd) for fit
    }
    #now map S^4 back to the original pairs:
    S2=rep(0,length(add))
    for (i in 1:length(L)){
	S2[L[[i]]]=S[i]^4
    }


    #no longer needed:
    # Calculate simple moderation term based only on LRV
    #mod <- z.df * z.s2 / propd@results$lrv
    
    mod <- z.df * S2 / propd@results$lrv #moderation term based on ratio mean-variance trend

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