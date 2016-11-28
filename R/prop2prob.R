#' Calculate Probability from Proportionality
#'
#' This experimental helper function calculates probability from
#'  proportionality. When supplying one \code{propr} object,
#'  \code{prop2prob} estimates the probability that
#'  each proportionality coefficient occurred by chance alone.
#'  When supplying two \code{propr} objects, \code{prop2prob}
#'  estimates the probability that each proportionality
#'  coefficient differs between the two objects.
#'
#' All calculations use formulae derived for the concordance
#'  correlation coefficient under the constraint that all means
#'  equal zero. We defend this constraint on the grounds that
#'  we can shift the mean of log-ratio transformed feature vectors
#'  without changing the proportionality coefficient, rho.
#'  See Zar's Biostatistical Analysis text (4ed, pg 407-10)
#'  for more information on the method used.
#'
#' When calculating differential proportionality, it is the
#'  responsibility of the user to ensure that the two groups
#'  have no overlapping samples. All p-values returned as
#'  twice the result of \code{\link{pnorm}}, thereby correcting
#'  for "two-tails". Please make sure to interpret p-values
#'  in the context of multiple testing! For more information,
#'  see \code{\link{p.adjust}}.
#'
#' @return A \code{data.table} of p-values.
#'
#' @param x A \code{propr} object.
#' @param y A \code{propr} object. Optional.
#'
#' @export
prop2prob <- function(x, y){

  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Uh oh! This display method depends on data.table! ",
         "Try running: install.packages('data.table')")
  }

  if(class(x) != "propr") stop("Uh oh! This function requires a 'propr' object for 'x'.")
  if(!x@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
  if(!missing(y)){
    if(class(y) != "propr") stop("Uh oh! This function requires a 'propr' object for 'y'.")
    if(!y@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
    if(!identical(colnames(x@logratio), colnames(y@logratio))){
      stop("Uh oh! Make sure both 'propr' objects have the same features.")
    }
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

  labels <- labRcpp(ncol(x@logratio))
  data.table::data.table(
    "Partner" = labels[[1]],
    "Pair" = labels[[2]],
    "Probability" = z
  )
}
