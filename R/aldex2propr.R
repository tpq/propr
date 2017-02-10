#' Import \code{ALDEx2} Object
#'
#' This method constructs a \code{propr} object from an
#'  \code{aldex.clr} object. See Details.
#'
#' The \code{ALDEx2} package has two exceptional features useful
#'  in proportionality analysis too. First, \code{ALDEx2} offers
#'  a number of additional log-ratio transformations, toggled
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
#' @param aldex.clr An \code{aldex.clr} object.
#' @param how A character string. The proportionality method
#'  used to build the \code{propr} object. For example,
#'  "perb" returns an estimation of rho while "phit" returns
#'  an estimation of phi.
#' @return Returns a \code{propr} object.
#'
#' @export
aldex2propr <- function(aldex.clr, how = "perb"){

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
  mc <- ALDEx2::getMonteCarloInstances(aldex.clr)
  k <- ALDEx2::numMCInstances(aldex.clr)
  logratio <- 0
  prop <- 0
  for(i in 1:k){

    cat(paste0(i, "..."))
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    logratio <- logratio + mci_lr
    prop.i <- do.call(how, list("lr" = mci_lr))
    prop <- prop + prop.i
  }
  cat("\n")

  # Clean up NAs in the case of denom = 1
  prop[is.na(prop)] <- ifelse(how == "lr2rho", 1, 0)

  propr <- new("propr")
  propr@counts <- t(as.matrix(aldex.clr@reads))
  propr@logratio <- logratio / k
  propr@matrix <- prop / k

  return(propr)
}
