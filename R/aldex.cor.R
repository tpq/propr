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
