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
#' @param how A character string. The proportionality method
#'  used to build the \code{propr} object. For example,
#'  "perb" returns an estimation of rho while "phit" returns
#'  an estimation of phi.
#' @param select See \code{\link{proportionality}}.
#' @return Returns a \code{propr} object.
#'
#' @export
aldex2propr <- function(aldex.clr, how = "perb", select){

  packageCheck("ALDEx2")

  if(class(aldex.clr) != "aldex.clr"){

    stop("This method expects an 'aldex.clr' object.")
  }

  if(how %in% c("perb", "rho", "lr2rho")){

    how <- "lr2rho"

  }else if(how %in% c("phit", "phi", "lr2phi")){

    how <- "lr2phi"

  }else if(how %in% c("phis", "phis", "phs", "lr2phs")){

    how <- "lr2phs"

  }else{

    stop("Provided 'how' not recognized.")
  }

  # Keep a running sum of propr instances
  counts <- t(as.matrix(aldex.clr@reads))
  mc <- ALDEx2::getMonteCarloInstances(aldex.clr)
  k <- ALDEx2::numMCInstances(aldex.clr)
  logratio <- 0
  prop <- 0
  for(i in 1:k){

    numTicks <- progress(i, k, numTicks)

    # Extract i-th Monte Carlo instance
    mci_lr <- t(sapply(mc, function(x) x[, i]))

    # Subset log-ratio transformed data
    if(!missing(select)){

      if(i == 1){

        # Make select boolean (it's OK if it's integer)
        if(is.character(select)) select <- match(select, colnames(mci_lr))
        if(any(is.na(select))) stop("Uh oh! Provided select reference not found in data.")
        counts <- counts[, select]
      }

      mci_lr <- mci_lr[, select]
    }

    # Add i-th log-ratio transformation to cumulative sum
    logratio <- logratio + mci_lr

    # Add i-th proportionality matrix to cumulative sum
    prop.i <- do.call(how, list("lr" = mci_lr))
    prop <- prop + prop.i
  }

  propr <- new("propr")
  propr@counts <- counts
  propr@logratio <- logratio / k
  propr@matrix <- prop / k

  return(propr)
}

#' Correlate CLR with a Continuous Measurement
#'
#' This function correlates each log-ratio transformed
#'  feature vector with a continuous numeric variable.
#'
#' @param lr A data.frame. A log-ratio transformed counts matrix.
#' @param conditions A numeric vector of a continuous variable.
#' @param ... Arguments passed to \code{cor.test}.
#'
#' @return Returns a data.frame of the correlation
#'  statistic (e.g., \code{r}) and p-value (\code{p})
#'  for each log-ratio transformed feature.
#'  FDR provided by \code{BH} column.
#'
#' @importFrom stats cor.test p.adjust
#' @export
lr2cor <- function(lr, conditions, ...){

  if(!is.numeric(conditions)){

    stop("Please provide a numeric 'conditions' argument.")
  }

  if(nrow(lr) != length(conditions)){

    stop("Incorrect length for 'conditions' argument.")
  }

  cors <- apply(lr, 2, function(x){
    cor.test(x, conditions, ...)
  })

  r <- sapply(cors, getElement, "statistic")
  p <- sapply(cors, getElement, "p.value")
  BH <- p.adjust(p, method = "BH")

  data.frame(r, p, BH,
             row.names = colnames(lr))
}

#' Correlate CLR with a Continuous Measurement
#'
#' This function uses the Monte Carlo instances from an
#'  \code{aldex.clr} object to correlate each log-ratio
#'  transformed feature vector with a continuous
#'  numeric variable. See \code{\link{lr2cor}}.
#'
#' @param clr An \code{aldex.clr} object.
#' @param conditions A numeric vector of a continuous variable.
#' @param ... Arguments passed to \code{cor.test}.
#'
#' @return Returns a data.frame of the average
#'  correlation statistic (e.g., \code{r}) and
#'  average p-value (\code{p}) for each feature
#'  across all Monte Carlo instances. Average
#'  FDR provided by \code{BH} column.
#'
#' @export
aldex.cor <- function(clr, conditions, ...){

  packageCheck("ALDEx2")

  # Keep a running sum of lr2cor instances
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for(i in 1:k){

    numTicks <- progress(i, k, numTicks)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2cor(mci_lr, conditions, ...)
  }

  r <- r / k
  data.frame("r" = r$r, "p" = r$p, "BH" = r$BH,
             row.names = colnames(mci_lr))
}
