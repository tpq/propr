#' Import \code{ALDEx2} Object
#'
#' This method constructs a \code{propr} object from an
#'  \code{aldex.clr} object. See Details.
#'
#' The \code{ALDEx2} package has two exceptional features useful
#'  in proportionality analysis too. First, \code{ALDEx2} offers
#'  a number of extra log-ratio transformations, toggled
#'  by the \code{denom} argument in \code{aldex.clr}. Second,
#'  \code{ALDEx2} estimates per-feature technical variation
#'  within each sample using Monte-Carlo instances drawn
#'  from the Dirichlet distribution.
#'
#' The \code{aldex2propr} function takes advantage of both
#'  of these features by constructing a \code{propr} object
#'  directly from an \code{aldex.clr} object. When interpreting
#'  the resultant \code{propr} object, keep in mind that
#'  \code{ALDEx2} adds 0.5 to all \code{@@counts} regardless
#'  of whether the counts contain any zeros. Otherwise,
#'  the \code{@@logratio} slot contains the log-ratio
#'  transformed counts as averaged across all Monte Carlo
#'  instances. Likewise, the \code{@@matrix} slot gets
#'  filled with the proportionality matrix as averaged
#'  across all Monte Carlo instances.
#'
#' The \code{select} argument subsets the feature matrix
#'  after log-ratio transformation but before calculating
#'  proportionality. This reduces the run-time and RAM
#'  overhead without impacting the final result. Removing
#'  lowly abundant features prior to log-ratio transformation
#'  could otherwise change the proportionality measure.
#'
#' @param aldex.clr An \code{aldex.clr} object.
#' @param how A character string. The proportionality metric
#'  used to build the \code{propr} object. Choose from
#'  "rho", "phi", or "phs".
#' @param select Optional. Use this to subset the final
#'  proportionality matrix without altering the result.
#'
#' @return Returns a \code{propr} object.
#'
#' @export
aldex2propr <- function(aldex.clr, how = "rho", select) {
  NVTX_PUSH("aldex2propr", 0)
  packageCheck("ALDEx2")

  if (!inherits(aldex.clr, "aldex.clr")) {
    NVTX_POP()  # aldex2propr
    stop("This method expects an 'aldex.clr' object.")
  }

  NVTX_PUSH("resolve_how", 0)
  if (how %in% c("perb", "rho", "lr2rho")) {
    how <- "lr2rho"
  } else if (how %in% c("phit", "phi", "lr2phi")) {
    how <- "lr2phi"
  } else if (how %in% c("phis", "phis", "phs", "lr2phs")) {
    how <- "lr2phs"
  } else {
    NVTX_POP()  # resolve_how
    NVTX_POP()  # aldex2propr
    stop("Provided 'how' not supported.")
  }
  NVTX_POP()  # resolve_how

  # Keep a running sum of propr instances
  NVTX_PUSH("init_counts_and_mc", 0)
  counts <- t(as.matrix(aldex.clr@reads))
  mc <- ALDEx2::getMonteCarloInstances(aldex.clr)
  k <- ALDEx2::numMCInstances(aldex.clr)
  logratio <- 0
  prop <- 0
  NVTX_POP()  # init_counts_and_mc

  for (i in 1:k) {
    NVTX_PUSH("iteration", 0)
    numTicks <- progress(i, k, numTicks)

    # Extract i-th Monte Carlo instance
    NVTX_PUSH("extract_MC_instance", 0)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    NVTX_POP()  # extract_MC_instance

    # Subset log-ratio transformed data
    if (!missing(select)) {
      NVTX_PUSH("subset_logratio", 0)
      if (i == 1) {
        # Make select boolean (it's OK if it's integer)
        if (is.character(select))
          select <- match(select, colnames(mci_lr))
        if (any(is.na(select))) {
          NVTX_POP()  # subset_logratio
          NVTX_POP()  # iteration
          NVTX_POP()  # aldex2propr
          stop("Uh oh! Provided select reference not found in data.")
        }
        counts <- counts[, select]
      }

      mci_lr <- mci_lr[, select]
      NVTX_POP()  # subset_logratio
    }

    # Add i-th log-ratio transformation to cumulative sum
    NVTX_PUSH("logratio_cumsum", 0)
    logratio <- logratio + mci_lr
    NVTX_POP()  # logratio_cumsum

    # Add i-th proportionality matrix to cumulative sum
    NVTX_PUSH("proportionality_cumsum", 0)
    prop.i <- do.call(how, list("lr" = mci_lr))
    prop <- prop + prop.i
    NVTX_POP()  # proportionality_cumsum

    NVTX_POP()  # iteration
  }

  NVTX_PUSH("build_propr_object", 0)
  propr <- new("propr")
  propr@counts <- as.data.frame(counts)
  propr@logratio <- as.data.frame(logratio) / k
  propr@matrix <- prop / k

  message("Alert: Using 'aldex2propr' is not compatible the @results table.")
  propr@results <- data.frame()

  message("Alert: Using 'aldex2propr' disables permutation testing.")
  propr@permutes <- list(NULL)
  NVTX_POP()  # build_propr_object

  NVTX_POP()  # aldex2propr
  return(propr)
}
