#' Update FDR by Permutation
#'
#' This function updates FDR for a set of cutoffs.
#'
#' This function wraps \code{updateCutoffs.propr} and
#'  \code{updateCutoffs.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param cutoff For \code{updateCutoffs}, a numeric vector.
#'  this argument provides the FDR cutoffs to test.
#' @param ncores An integer. The number of parallel cores to use.
#' @return A \code{propr} or \code{propd} object with the FDR slot updated.
#' @export
updateCutoffs <-
  function(object,
           cutoff = seq(.05, .95, .3),
           ncores = 1) {
    if (inherits(object, "propr")) {
      if (ncores == 1) {
        message("Alert: Try parallelizing updateCutoffs with ncores > 1.")
      }

      updateCutoffs.propr(object, cutoff, ncores)

    } else if (inherits(object, "propd")) {
      if (ncores > 1) {
        message("Alert: Parallel updateCutoffs not yet supported.")
      }

      updateCutoffs.propd(object, cutoff)

    } else{
      stop("Provided 'object' not recognized.")
    }
  }

#' @rdname updateCutoffs
#' @section Methods:
#' \code{updateCutoffs.propr:}
#'  Use the \code{propr} object to permute proportionality
#'  across a number of cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#' @export
updateCutoffs.propr <-
  function(object, cutoff, ncores) {

    metric_is_up <- function(metric) {
      metrics <- c("rho", "cor", "pcor", "pcor.shrink", "pcor.bshrink")
      return(metric %in% metrics)
    }

    getFdrRandcounts <- function(ct.k) {
      pr.k <- suppressMessages(propr::propr(
        ct.k,
        object@metric,
        ivar = object@ivar,
        alpha = object@alpha,
        p = 0
      ))

      # Vector of propr scores for each pair of taxa.
      pkt <- pr.k@results$propr

      # Find number of permuted theta less than cutoff
      sapply(FDR$cutoff, function(cut) {
        if (metric_is_up(object@metric)) {
          count_greater_than(pkt, cut)
        } else{
          count_less_than(pkt, cut)
        }
      })
    }

    if (object@metric == "rho") {
      message("Alert: Estimating FDR for largely positive proportional pairs only.")
    }

    if (object@metric == "phi") {
      warning("We recommend using the symmetric phi 'phs' for FDR permutation.")
    }

    if (identical(object@permutes, list(NULL)))
      stop("Permutation testing is disabled.")

    # Let NA cutoff skip function
    if (identical(cutoff, NA))
      return(object)

    # Set up FDR cutoff table
    FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
    FDR$cutoff <- cutoff
    p <- length(object@permutes)

    if (ncores > 1) {
      packageCheck("parallel")

      # Set up the cluster and require propr
      cl <- parallel::makeCluster(ncores)
      # parallel::clusterEvalQ(cl, requireNamespace(propr, quietly = TRUE))

      # Each element of this list will be a vector whose elements
      # are the count of theta values less than the cutoff.
      randcounts <- parallel::parLapply(cl = cl,
                                        X = object@permutes,
                                        fun = getFdrRandcounts)

      # Sum across cutoff values
      FDR$randcounts <- apply(as.data.frame(randcounts), 1, sum)

      # Explicitly stop the cluster.
      parallel::stopCluster(cl)

    } else{
      # Calculate propr for each permutation -- NOTE: `select` and `subset` disable permutation testing
      for (k in 1:p) {
        numTicks <- progress(k, p, numTicks)

        # Calculate propr exactly based on @metric, @ivar, and @alpha
        ct.k <- object@permutes[[k]]
        pr.k <- suppressMessages(propr(
          ct.k,
          object@metric,
          ivar = object@ivar,
          alpha = object@alpha,
          p = 0
        ))
        pkt <- pr.k@results$propr

        # Find number of permuted theta less than cutoff
        for (cut in 1:nrow(FDR)) {
          # randcounts as cumsum

          # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
          if (metric_is_up(object@metric)) {
            FDR[cut, "randcounts"] <-
              FDR[cut, "randcounts"] + count_greater_than(pkt, FDR[cut, "cutoff"])
          } else{
            # phi & phs
            FDR[cut, "randcounts"] <-
              FDR[cut, "randcounts"] + count_less_than(pkt, FDR[cut, "cutoff"])
          }
        }
      }
    }

    # Calculate FDR based on real and permuted tallys
    FDR$randcounts <- FDR$randcounts / p # randcounts as mean

    for (cut in 1:nrow(FDR)) {
      # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
      if (metric_is_up(object@metric)) {
        FDR[cut, "truecounts"] <-
          sum(object@results$propr > FDR[cut, "cutoff"])
      } else{
        # phi & phs
        FDR[cut, "truecounts"] <-
          sum(object@results$propr < FDR[cut, "cutoff"])
      }

      FDR[cut, "FDR"] <-
        FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
    }

    # Initialize @fdr
    object@fdr <- FDR

    return(object)
  }

#' @rdname updateCutoffs
#' @section Methods:
#' \code{updateCutoffs.propd:}
#'  Use the \code{propd} object to permute theta across a
#'  number of theta cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#' @export
updateCutoffs.propd <-
  function(object, cutoff = seq(.05, .95, .3)) {
    if (identical(object@permutes, data.frame()))
      stop("Permutation testing is disabled.")

    # Let NA cutoff skip function
    if (identical(cutoff, NA))
      return(object)

    # Set up FDR cutoff table
    FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
    FDR$cutoff <- cutoff
    p <- ncol(object@permutes)
    lrv <- object@results$lrv

    # Use calculateTheta to permute active theta
    for (k in 1:p) {
      numTicks <- progress(k, p, numTicks)

      # Tally k-th thetas that fall below each cutoff
      shuffle <- object@permutes[, k]

      if (object@active == "theta_mod") {
        # Calculate theta_mod with updateF (using i-th permuted object)
        if (is.na(object@Fivar))
          stop("Please re-run 'updateF' with 'moderation = TRUE'.")
        propdi <- suppressMessages(
          propd(
            object@counts[shuffle,],
            group = object@group,
            alpha = object@alpha,
            p = 0,
            weighted = object@weighted
          )
        )
        propdi <-
          suppressMessages(updateF(propdi, moderated = TRUE, ivar = object@Fivar))
        pkt <- propdi@results$theta_mod

      } else{
        # Calculate all other thetas directly (using calculateTheta)
        pkt <- suppressMessages(
          calculate_theta(
            object@counts[shuffle,],
            object@group,
            object@alpha,
            lrv,
            only = object@active,
            weighted = object@weighted
          )
        )
      }

      # Find number of permuted theta less than cutoff
      for (cut in 1:nrow(FDR)) {
        # randcounts as cumsum
        FDR[cut, "randcounts"] <-
          FDR[cut, "randcounts"] + sum(pkt < FDR[cut, "cutoff"])
      }
    }

    # Calculate FDR based on real and permuted tallys
    FDR$randcounts <- FDR$randcounts / p # randcounts as mean
    for (cut in 1:nrow(FDR)) {
      FDR[cut, "truecounts"] <-
        sum(object@results$theta < FDR[cut, "cutoff"])
      FDR[cut, "FDR"] <-
        FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
    }

    # Initialize @fdr
    object@fdr <- FDR

    return(object)
  }
