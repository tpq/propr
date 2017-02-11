#' Correlate CLR with a Continuous Measurement
#'
#' This function correlates each log-ratio transformed
#'  feature vector with a continuous numeric variable.
#'
#' @param lr A data.frame. A log-ratio transformed counts matrix.
#' @param conditions A numeric vector of a continuous variable.
#'
#' @return Returns a data.frame of the Pearson's correlation
#'  coefficient (\code{r}), Fisher's z-transformation (\code{z}),
#'  and p-value (\code{p}) for each feature.
#'
#' @export
lr2cor <- function(lr, conditions){

  if(!is.numeric(conditions)){

    stop("Please provide a numeric 'conditions' argument.")
  }

  if(nrow(lr) != length(conditions)){

    stop("Incorrect length for 'conditions' argument.")
  }

  r <- apply(lr, 2, cor, conditions)
  z <- atanh(r)
  sd <- 1/sqrt(length(conditions) - 3)
  p <- pnorm(abs(z), sd = sd, lower.tail = FALSE)

  data.frame(r, z, p)
}

#' Correlate CLR with a Continuous Measurement
#'
#' This function uses the Monte Carlo instances from an
#'  \code{aldex.clr} object to correlate each log-ratio
#'  transformed feature vector with a continuous
#'  numeric variable.
#'
#' @param clr An \code{aldex.clr} object.
#' @param conditions A numeric vector of a continuous variable.
#'
#' @return Returns a data.frame of the average
#'  Pearson's correlation coefficient (\code{r}),
#'  a Fisher's z-transformation (\code{z}) of this
#'  average, and the associated p-value (\code{p})
#'  for each feature across all Monte Carlo
#'  instances.
#'
#' @export
aldex.cor <- function(clr, conditions){

  packageCheck("ALDEx2")

  # Keep a running sum of lr2cor instances
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for(i in 1:k){

    cat(paste0(i, "..."))
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2cor(mci_lr, conditions)$r
  }
  cat("\n")

  r <- r / k
  z <- atanh(r)
  sd <- 1/sqrt(length(conditions) - 3)
  p <- pnorm(abs(z), sd = sd, lower.tail = FALSE)

  data.frame(r, z, p,
             row.names = colnames(mci_lr))
}
