#' Calculate Probability from Proportionality
#'
#' This experimental helper function calculates probability from
#'  proportionality. When supplying one \code{propr} object,
#'  \code{prop2prob} estimates the probability that
#'  each proportionality coefficient occurred by chance alone.
#'  When supplying two \code{propr} objects, \code{prop2prob}
#'  estimates the probability that each proportionality
#'  coefficient differs between the two groups, called
#'  differential proportionality.
#'
#' All calculations use formula derived for the concordance
#'  correlation coefficient under the constraint that all means
#'  equal zero. We defend this constraint on the grounds that
#'  we can shift the mean of log-ratio transformed feature vectors
#'  without changing the proportionality coefficient, rho.
#'  See Zar's Biostatistical Analysis text (4ed, pg 407)
#'  for more information.
#'
#' When calculating differential proportionality, it is the
#'  responsibility of the user to ensure that the two groups
#'  have no overlapping samples. All p-values calculated as
#'  twice the result of \code{\link{pnorm}}, thereby correcting
#'  for "two-tails". Please make sure to interpret p-values
#'  in the context of multiple testing! For more information,
#'  see \code{\link{p.adjust}}.
#'
#' @param x A \code{propr} object.
#' @param y A \code{propr} object. Optional.
#'
#' @export
prop2prob <- function(x, y){

  if(class(x) != "propr") stop("Uh oh! This function requires a 'propr' object for 'x'.")
  if(!x@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
  N.x <- ncol(x@logratio)
  if(!missing(y)){
    if(class(y) != "propr") stop("Uh oh! This function requires a 'propr' object for 'y'.")
    if(!y@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
    if(N.x != ncol(y@logratio)) stop("Uh oh! All objects must have the same number of features.")
  }

  X <- linRcpp(x@matrix, x@logratio)
  z <- lltRcpp(X)
  var <- urtRcpp(X)
  rm(X); gc()

  if(!missing(y)){

    Y <- linRcpp(y@matrix, y@logratio)
    z <- z - lltRcpp(Y)
    var <- var + urtRcpp(Y)
    rm(Y); gc()
  }

  # Calculate normal deviate
  z <- abs(z / sqrt(var))
  rm(var); gc()

  # Calculate probability
  z <- pnorm(z, lower.tail = FALSE) * 2

  labels <- labRcpp(N.x)
  data.frame(
    "Partner" = labels[[1]],
    "Pair" = labels[[2]],
    "Probability" = z
  )
}
