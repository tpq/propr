#' Update FDR by Permutation
#'
#' This function updates the FDR for each cutoff. By default, the set of cutoffs are determined
#' based on the quantile of the data, so that the cutoffs are evenly spaced across the data.
#' The FDR is calculated as the ratio between the number of permuted values beyond the cutoff
#' and the number of true values beyond the the cutoff. 
#' When tails is set to 'right', the FDR is calculated only on the positive side of the data.
#' When tails is set to 'both', the FDR is calculated on both sides of the data.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param number_of_cutoffs An integer. The number of cutoffs to test. Given this number, 
#' the cutoffs will be determined based on the quantile of the data. In this way, the 
#' cutoffs will be evenly spaced across the data.
#' @param custom_cutoffs A numeric vector. When provided, this vector is used as the set of 
#' cutoffs to test, and 'number_of_cutoffs' is ignored.
#' @param tails 'right' or 'both'. 'right' is for one-sided on the right. 'both' for
#' symmetric two-sided test. This is only relevant for \code{propr} objects, as 
#' \code{propd} objects are always one-sided and only have positive values. Default 
#' is 'right'.
#' @param ncores An integer. The number of parallel cores to use.
#' @return A \code{propr} or \code{propd} object with the FDR slot updated.
#' 
#' @export
updateCutoffs <-
  function(object,
           number_of_cutoffs = 100,
           custom_cutoffs = NA,
           tails = 'right',
           ncores = 1) {

    NVTX_PUSH("updateCutoffs", 0)

    if (inherits(object, "propr")) {
      res <- updateCutoffs.propr(object, number_of_cutoffs, custom_cutoffs, tails, ncores)
      NVTX_POP()
      return(res)

    } else if (inherits(object, "propd")) {
      res <- updateCutoffs.propd(object, number_of_cutoffs, custom_cutoffs, ncores)
      NVTX_POP()
      return(res)

    } else {
      NVTX_POP()
      stop("Provided 'object' not recognized.")
    }

  }

#' @rdname updateCutoffs
#' @section Methods:
#' \code{updateCutoffs.propr:}
#'  Use the \code{propr} object to permute correlation-like metrics
#'  (ie. rho, phi, phs, cor, pcor, pcor.shrink, pcor.bshrink),
#'  across a number of cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#' @export
updateCutoffs.propr <-
  function(object, 
           number_of_cutoffs = 100, 
           custom_cutoffs = NA, 
           tails = 'right', 
           ncores = 1) {

    if (identical(object@permutes, list(NULL))) {
      stop("Permutation testing is disabled.")
    }

    NVTX_PUSH("updateCutoffs.propr", 0)

    if (object@metric == "phi") {
      warning("We recommend using the symmetric phi 'phs' for FDR permutation.")
    }
    object@tails <- tails

    # get cutoffs
    NVTX_PUSH("get_cutoffs.propr", 0)
    if (length(custom_cutoffs) == 1 && is.na(custom_cutoffs)) {
      vals <- object@results$propr
      if (tails == 'right') {
        vals <- vals[vals >= 0]
      } else if (tails == 'both') {
        vals <- abs(vals)
      }
      cutoffs <- as.numeric(quantile(vals, seq(0, 1, length.out = number_of_cutoffs)))
    } else {
      cutoffs <- custom_cutoffs
    }
    NVTX_POP()  # get_cutoffs.propr

    # Set up FDR cutoff table
    NVTX_PUSH("setup_FDR_table.propr", 0)
    FDR <- as.data.frame(matrix(0, nrow = length(cutoffs), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
    FDR$cutoff <- cutoffs
    NVTX_POP()  # setup_FDR_table.propr

    # count the permuted values greater or less than each cutoff
    NVTX_PUSH("randcounts.propr", 0)
    if (ncores > 1) {
      FDR$randcounts <- getFdrRandcounts.propr.parallel(object, cutoffs, ncores)
    } else {
      FDR$randcounts <- getFdrRandcounts.propr.run(object, cutoffs)
    }
    NVTX_POP()  # randcounts.propr

    # count actual values greater or less than each cutoff
    NVTX_PUSH("truecounts.propr", 0)
    vals <- object@results$propr
    if (tails == 'both') vals <- abs(vals)
    FDR$truecounts <- sapply(FDR$cutoff, function(cutoff) {
      countValuesBeyondThreshold(vals, cutoff, object@direct)
    })
    NVTX_POP()  # truecounts.propr

    # calculate FDR
    NVTX_PUSH("calc_FDR.propr", 0)
    FDR$FDR <- FDR$randcounts / FDR$truecounts
    object@fdr <- FDR
    NVTX_POP()  # calc_FDR.propr

    NVTX_POP()  # updateCutoffs.propr
    return(object)
  }

#' This function counts the permuted values greater or less than each cutoff,
#' using parallel processing, for a propr object
getFdrRandcounts.propr.parallel <- 
  function(object, cutoffs, ncores) {

    NVTX_PUSH("getFdrRandcounts.propr.parallel", 0)

    # Set up the cluster
    NVTX_PUSH("setup_cluster.propr", 0)
    packageCheck("parallel")
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterEvalQ(cl, requireNamespace(propr, quietly = TRUE))
    NVTX_POP()  # setup_cluster.propr

    # define function to parallelize
    getFdrRandcounts <- function(ct.k) {
      # calculate permuted propr
      pr.k <- suppressMessages(propr::propr(
        ct.k,
        object@metric,
        ivar = object@ivar,
        alpha = object@alpha,
        p = 0
      ))
      # Vector of propr scores for each pair of taxa.
      pkt <- pr.k@results$propr
      if (object@tails == 'both') pkt <- abs(pkt)
      # Find number of permuted theta more or less than cutoff
      sapply(cutoffs, function(cutoff) countValuesBeyondThreshold(pkt, cutoff, object@direct))
    }

    # Each element of this list will be a vector whose elements
    # are the count of theta values less than the cutoff.
    NVTX_PUSH("parLapply.propr", 0)
    randcounts <- parallel::parLapply(cl = cl,
                                      X = object@permutes,
                                      fun = getFdrRandcounts)
    NVTX_POP()  # parLapply.propr

    # get the average randcounts across all permutations
    NVTX_PUSH("aggregate_randcounts.propr", 0)
    randcounts <- apply(as.data.frame(randcounts), 1, sum)
    randcounts <- randcounts / length(object@permutes)
    NVTX_POP()  # aggregate_randcounts.propr

    # Explicitly stop the cluster.
    NVTX_PUSH("stop_cluster.propr", 0)
    parallel::stopCluster(cl)
    NVTX_POP()  # stop_cluster.propr

    NVTX_POP()  # getFdrRandcounts.propr.parallel
    return(randcounts)
  }

#' This function counts the permuted values greater or less than each cutoff,
#' using a single core, for a propr object.
getFdrRandcounts.propr.run <-
  function(object, cutoffs) {

    NVTX_PUSH("getFdrRandcounts.propr.run", 0)

    # create empty randcounts
    randcounts <- rep(0, length(cutoffs))

    # Calculate propr for each permutation -- NOTE: `select` and `subset` disable permutation testing
    p <- length(object@permutes)

    NVTX_PUSH("permutation_loop.propr", 0)
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
      if (object@tails == 'both') pkt <- abs(pkt)

      # calculate the cumulative (across permutations) number of permuted values more or less than cutoff
      for (cut in 1:length(cutoffs)) {
        randcounts[cut] <- randcounts[cut] + countValuesBeyondThreshold(pkt, cutoffs[cut], object@direct)
      }
    }
    NVTX_POP()  # permutation_loop.propr

    randcounts <- randcounts / p  # averaged across permutations

    NVTX_POP()  # getFdrRandcounts.propr.run
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
  function(object, number_of_cutoffs = 100, custom_cutoffs = NA, ncores = 1) {

    if (identical(object@permutes, data.frame()))
      stop("Permutation testing is disabled.")

    NVTX_PUSH("updateCutoffs.propd", 0)

    # get cutoffs
    NVTX_PUSH("get_cutoffs.propd", 0)
    if (length(custom_cutoffs) == 1 && is.na(custom_cutoffs)) {
      cutoffs <- as.numeric(quantile(object@results$theta, seq(0, 1, length.out = number_of_cutoffs)))
    } else {
      cutoffs <- custom_cutoffs
    }
    NVTX_POP()  # get_cutoffs.propd

    # Set up FDR cutoff table
    NVTX_PUSH("setup_FDR_table.propd", 0)
    FDR <- as.data.frame(matrix(0, nrow = length(cutoffs), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
    FDR$cutoff <- cutoffs
    NVTX_POP()  # setup_FDR_table.propd

    # Count the permuted values greater or less than each cutoff
    NVTX_PUSH("randcounts.propd", 0)
    if (ncores > 1) {
      FDR$randcounts <- getFdrRandcounts.propd.parallel(object, cutoffs, ncores)
    } else {
      FDR$randcounts <- getFdrRandcounts.propd.run(object, cutoffs)
    }
    NVTX_POP()  # randcounts.propd

    # count actual values greater or less than each cutoff
    NVTX_PUSH("truecounts.propd", 0)
    FDR$truecounts <- sapply(1:nrow(FDR), function(cut) {
      countValuesBeyondThreshold(object@results$theta, FDR[cut, "cutoff"], direct = FALSE)
    })
    NVTX_POP()  # truecounts.propd

    # Calculate FDR
    NVTX_PUSH("calc_FDR.propd", 0)
    FDR$FDR <- FDR$randcounts / FDR$truecounts
    object@fdr <- FDR
    NVTX_POP()  # calc_FDR.propd

    NVTX_POP()  # updateCutoffs.propd
    return(object)
  }

#' This function counts the permuted values greater or less than each cutoff,
#' using parallel processing, for a propd object.
getFdrRandcounts.propd.parallel <- 
  function(object, cutoffs, ncores) {

    NVTX_PUSH("getFdrRandcounts.propd.parallel", 0)

    # Set up the cluster
    NVTX_PUSH("setup_cluster.propd", 0)
    packageCheck("parallel")
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterEvalQ(cl, requireNamespace(propr, quietly = TRUE))
    NVTX_POP()  # setup_cluster.propd

    # define functions to parallelize
    getFdrRandcountsMod <- function(k) {
      if (is.na(object@Fivar)) stop("Please re-run 'updateF' with 'moderation = TRUE'.")
      shuffle <- object@permutes[, k]
      propdi <- suppressMessages(
        propd(
          object@counts[shuffle,],
          group = object@group,
          alpha = object@alpha,
          p = 0,
          weighted = object@weighted,
          shrink = object@shrink
        )
      )
      propdi <- suppressMessages(updateF(propdi, moderated = TRUE, ivar = object@Fivar))
      pkt <- propdi@results$theta_mod
      sapply(cutoffs, function(cutoff) countValuesBeyondThreshold(pkt, cutoff, direct = FALSE))
    }
    getFdrRandcounts <- function(k) {
      shuffle <- object@permutes[, k]
      pkt <- suppressMessages(
        calculate_theta(
          object@counts[shuffle,],
          object@group,
          object@alpha,
          object@results$lrv,
          only = object@active,
          weighted = object@weighted,
          shrink = object@shrink
        )
      )
      sapply(cutoffs, function(cutoff) countValuesBeyondThreshold(pkt, cutoff, direct = FALSE))
    }

    # Each element of this list will be a vector whose elements
    # are the count of theta values less than the cutoff.
    NVTX_PUSH("parLapply.propd", 0)
    if (object@active == "theta_mod") {
      func <- getFdrRandcountsMod
    } else {
      func <- getFdrRandcounts
    }
    randcounts <- parallel::parLapply(cl = cl,
                                      X = 1:ncol(object@permutes),
                                      fun = func)
    NVTX_POP()  # parLapply.propd

    # get the average randcounts across all permutations
    NVTX_PUSH("aggregate_randcounts.propd", 0)
    randcounts <- apply(as.data.frame(randcounts), 1, sum)
    randcounts <- randcounts / length(object@permutes)
    NVTX_POP()  # aggregate_randcounts.propd

    # Explicitly stop the cluster.
    NVTX_PUSH("stop_cluster.propd", 0)
    parallel::stopCluster(cl)
    NVTX_POP()  # stop_cluster.propd

    NVTX_POP()  # getFdrRandcounts.propd.parallel
    return(randcounts)
  }

#' This function counts the permuted values greater or less than each cutoff,
#' using a single core, for a propd object.
getFdrRandcounts.propd.run <- 
  function(object, cutoffs) {

    NVTX_PUSH("getFdrRandcounts.propd.run", 0)

    # create empty randcounts
    randcounts <- rep(0, length(cutoffs))

    # use calculateTheta to permute active theta
    p <- ncol(object@permutes)

    NVTX_PUSH("permutation_loop.propd", 0)
    for (k in 1:p) {
      numTicks <- progress(k, p, numTicks)

      # calculate permuted theta values
      if (object@active == "theta_mod") {
        pkt <- suppressMessages(getPermutedThetaMod(object, k))
      } else {
        pkt <- suppressMessages(getPermutedTheta(object, k))
      }

      # calculate the cumulative (across permutations) number of permuted values more or less than cutoff
      for (cut in 1:length(cutoffs)) {
        randcounts[cut] <- randcounts[cut] + countValuesBeyondThreshold(pkt, cutoffs[cut], direct = FALSE)
      }
    }
    NVTX_POP()  # permutation_loop.propd

    randcounts <- randcounts / p  # averaged across permutations

    NVTX_POP()  # getFdrRandcounts.propd.run
    return(randcounts)
}

#' Get the theta mod values for a given permutation
getPermutedThetaMod <- 
  function(object, k) {

    if (is.na(object@Fivar)) stop("Please re-run 'updateF' with 'moderation = TRUE'.")

    NVTX_PUSH("getPermutedThetaMod", 0)

    # Tally k-th thetas that fall below each cutoff
    shuffle <- object@permutes[, k]

    # Calculate theta_mod with updateF (using k-th permuted object)
    propdi <- suppressMessages(
      propd(
        object@counts[shuffle,],
        group = object@group,
        alpha = object@alpha,
        p = 0,
        weighted = object@weighted,
        shrink = object@shrink
      )
    )
    propdi <- suppressMessages(updateF(propdi, moderated = TRUE, ivar = object@Fivar))

    res <- propdi@results$theta_mod
    NVTX_POP()  # getPermutedThetaMod
    return(res)
}

#' Get the theta values for a given permutation
getPermutedTheta <- 
  function(object, k) {

    NVTX_PUSH("getPermutedTheta", 0)

    # Tally k-th thetas that fall below each cutoff
    shuffle <- object@permutes[, k]

    # Calculate all other thetas directly (using calculateTheta)
    pkt <- suppressMessages(
      calculate_theta(
        object@counts[shuffle,],
        object@group,
        object@alpha,
        object@results$lrv,
        only = object@active,
        weighted = object@weighted,
        shrink = object@shrink
      )
    )

    NVTX_POP()  # getPermutedTheta
    return(pkt)
  }

#' Count Values Greater or Less Than a Threshold
#'
#' This function counts the number of values greater or less than a threshold.
#' The direction depends on if a direct or inverse relationship is asked,
#' as well as the sign of the threshold.
#'
#' @param values A numeric vector.
#' @param cutoff A numeric value.
#' @param direct A logical value. If \code{TRUE}, direct relationship is considered.
#' @return The number of values greater or less than the threshold.
countValuesBeyondThreshold <- function(values, cutoff, direct) {
  NVTX_PUSH("countValuesBeyondThreshold", 0)

  if (is.na(cutoff)) {
    NVTX_POP()
    return(NA)
  }

  if (direct && cutoff >= 0) {
    res <- count_greater_than(values, cutoff)
  } else {
    res <- count_less_than(values, cutoff)
  }

  NVTX_POP()
  return(res)
}
