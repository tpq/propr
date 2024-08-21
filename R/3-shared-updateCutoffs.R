#' Update FDR by Permutation
#'
#' This function updates FDR for a set of cutoffs.
#'
#' This function wraps \code{updateCutoffs.propr} and
#'  \code{updateCutoffs.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param number_of_cutoffs An integer. The number of cutoffs to test.
#' Given this number, the cutoffs will be determined based on the quantile of the data.
#' In this way, the cutoffs will be evenly spaced across the data.
#' @param custom_cutoffs A numeric vector. When provided, this vector is used
#' as the FDR cutoffs to test, and number_of_cutoffs is ignored.
#' @param ncores An integer. The number of parallel cores to use.
#' @return A \code{propr} or \code{propd} object with the FDR slot updated.
#' @export
updateCutoffs <-
  function(object,
           number_of_cutoffs = 100,
           custom_cutoffs = NULL,
           ncores = 1) {

    get_cutoffs <- function(values, number_of_cutoffs, custom_cutoffs) {
      if (!is.null(custom_cutoffs)) {
        return(custom_cutoffs)
      }else{
        return(as.numeric(quantile(values, probs = seq(0, 1, length.out = number_of_cutoffs))))
      }
    }

    if (inherits(object, "propr")) {
      values <- object@matrix[lower.tri(object@matrix)]
      cutoffs <- get_cutoffs(values, number_of_cutoffs, custom_cutoffs)
      updateCutoffs.propr(object, cutoffs, ncores)

    } else if (inherits(object, "propd")) {
      values <- object@results$theta
      cutoffs <- get_cutoffs(values, number_of_cutoffs, custom_cutoffs)
      updateCutoffs.propd(object, cutoffs)

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

    if (identical(object@permutes, list(NULL))) {
      stop("Permutation testing is disabled.")
    }
    if (object@metric == "rho") {
      message("Alert: Estimating FDR for largely positive proportional pairs only.")
    }
    if (object@metric == "phi") {
      warning("We recommend using the symmetric phi 'phs' for FDR permutation.")
    }

    # Set up FDR cutoff table
    FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
    FDR$cutoff <- cutoff

    # define the counting functions (greater or less than a threshold)
    # countFunc for counting positive values, and countFunNegative for negative values
    countFunc <- if (object@direct) count_greater_than else count_less_than
    countFunNegative <- if (object@direct) count_less_than else count_greater_than

    # count the permuted values greater or less than each cutoff
    if (ncores > 1) {
      FDR$randcounts <- updateCutoffs.propr.parallel(object, FDR$cutoff, ncores, countFunc, countFunNegative)
    } else{
      FDR$randcounts <- updateCutoffs.propr.run(object, FDR$cutoff, countFunc, countFunNegative)
    }

    # count actual values greater or less than each cutoff
    for (cut in 1:nrow(FDR)){
      if (FDR[cut, "cutoff"] > 0) currentFunc = countFunc else currentFunc = countFunNegative
      FDR[cut, "truecounts"] <- currentFunc(object@results$propr, FDR[cut, "cutoff"])
    }

    # calculate FDR
    FDR$FDR <- FDR$randcounts / FDR$truecounts

    # Initialize @fdr
    object@fdr <- FDR

    return(object)
  }

# run updateCutoffs.propr in parallel
updateCutoffs.propr.parallel <- 
  function(object, cutoffs, ncores, countFunc, countFunNegative) {

    # Set up the cluster
    packageCheck("parallel")
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterEvalQ(cl, requireNamespace(propr, quietly = TRUE))

    # define the function to parallelize
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
      # Find number of permuted theta more or less than cutoff
      sapply(cutoffs, function(cut) if (cut > 0) countFunc(pkt, cut) else countFunNegative(pkt, cut))
    }

    # Each element of this list will be a vector whose elements
    # are the count of theta values less than the cutoff.
    randcounts <- parallel::parLapply(cl = cl,
                                      X = object@permutes,
                                      fun = getFdrRandcounts)

    # Sum across cutoff values
    randcounts <- apply(as.data.frame(randcounts), 1, sum)
    randcounts <- randcounts / length(object@permutes)

    # Explicitly stop the cluster.
    parallel::stopCluster(cl)

    return(randcounts)
  }

# run updateCutoffs.propr not in parallel
updateCutoffs.propr.run <-
  function(object, cutoffs, countFunc, countFunNegative) {

    # Calculate propr for each permutation -- NOTE: `select` and `subset` disable permutation testing
    p <- length(object@permutes)
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

      # Find number of permuted theta more or less than cutoff
      randcounts <- rep(0, length(cutoffs))
      for (cut in 1:length(cutoffs)){
        if (cutoffs[cut] > 0) currentFunc = countFunc else currentFunc = countFunNegative
        randcounts[cut] <- randcounts[cut] + currentFunc(pkt, cutoffs[cut])
      }
    }

    randcounts <- randcounts / p
    return(randcounts)
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
  function(object, cutoff) {
    if (identical(object@permutes, data.frame()))
      stop("Permutation testing is disabled.")

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
      FDR$randcounts <- sapply(
        1:nrow(FDR), 
        function(cut) FDR$randcounts[cut] + count_less_than(pkt, FDR[cut, "cutoff"]),
        simplify = TRUE
      )
    }

    # Calculate FDR based on real and permuted tallys
    FDR$randcounts <- FDR$randcounts / p # randcounts as mean
    FDR$truecounts <- sapply(
      1:nrow(FDR), 
      function(cut) count_less_than(object@results$theta, FDR[cut, "cutoff"]), 
      simplify = TRUE
    )
    FDR$FDR <- FDR$randcounts / FDR$truecounts

    # Initialize @fdr
    object@fdr <- FDR

    return(object)
  }
