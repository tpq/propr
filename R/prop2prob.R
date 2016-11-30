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
#'  without changing the proportionality coefficient, rho, or
#'  Pearson's correlation coefficient, r. We refer the reader to
#'  Zar's Biostatistical Analysis text (4ed, pg 407-10) for
#'  more information on the method used.
#'
#' When calculating differential proportionality, it is the
#'  responsibility of the user to ensure that the two groups
#'  have no overlapping samples. All p-values returned as
#'  twice the result of \code{\link{pnorm}}, thereby correcting
#'  for "two-tails". Please make sure to interpret p-values
#'  in the context of multiple testing! For more information,
#'  see \code{\link{p.adjust}}.
#'
#' @param x A \code{propr} object.
#' @param y A \code{propr} object. Optional.
#' @param method A character string. Selects method used to
#'  adjust p-values for multiple comparisons. Argument
#'  passed to \code{\link{p.adjust}}. Defaults to the
#'  more conservative Bonferroni correction.
#' @inheritParams bucket
#'
#' @return A \code{data.table} of p-values.
#'
#' @seealso \code{\link{propr}}, \code{\link{abstract}}
#'
#' @examples
#' library(propr)
#' data(mail)
#' rho <- perb(mail)
#' prop2prob(rho)
#' @importFrom stats pnorm p.adjust
#' @export
prop2prob <- function(x, y, method = "bonferroni", prompt = TRUE){

  if(object.size(x) > 10^9 & prompt){

    message("Uh oh! A large proportionality matrix was detected (> 1 GB).\n",
            "This operation requires about 5 times the size of 'x' in additional RAM.\n",
            "Are you sure you want to calculate ",
            ncol(rho@matrix) * (ncol(rho@matrix)-1) * 1/2, " p-values?\n",
            "0: Nevermind\n1: Proceed\n2: Hmm...")
    response <- readline(prompt = "Which do you choose? ")
    if(!response == 1) stop("prop2prob method aborted.")
  }

  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Uh oh! This display method depends on data.table! ",
         "Try running: install.packages('data.table')")
  }

  differentialCheck(x, y, forceBoth = FALSE)

  X <- linRcpp(x@matrix[], x@logratio[])
  z <- lltRcpp(X)
  var <- urtRcpp(X)
  rm(X); gc()

  if(!missing(y)){

    Y <- linRcpp(y@matrix[], y@logratio[])
    z <- z - lltRcpp(Y)
    var <- var + urtRcpp(Y)
    rm(Y); gc()
  }

  # Calculate normal deviate
  z <- abs(z / sqrt(var))
  rm(var); gc()

  # Calculate probability
  z <- pnorm(z, lower.tail = FALSE) * 2
  a <- p.adjust(z, method = method)

  labels <- labRcpp(ncol(x@logratio))
  dt <- data.table::data.table(
    "Partner" = labels[[1]],
    "Pair" = labels[[2]],
    "Probability" = z,
    "Adjusted" = a,
    key = "Probability"
  )

  return(dt)
}

#' Abstract Two propr Objects
#'
#' This function abstracts a new \code{propr} object from
#'  two existing \code{propr} objects. The two \code{propr}
#'  objects should not have any overlapping samples. Typically,
#'  the two objects represent different experimental groups.
#'  The resultant abstracted object inherits all plot functions
#'  available for the original \code{propr} objects.
#'
#' The abstracted \code{propr} object has the following properties:
#'  The \code{@@counts} and \code{@@logratio} slots contain a
#'  join of the original slots via \code{rbind}. Meanwhile,
#'  the \code{@@matrix} slot contains a difference matrix defined as
#'  \code{tanh(atanh(x@matrix) - atanh(y@matrix))}. Viewing this difference
#'  matrix with \code{\link{dendrogram}} may help summarize
#'  the results of \code{\link{prop2prob}}. Note that this
#'  difference matrix now also informs co-cluster assignment for
#'  \code{\link{bucket}}, \code{\link{prism}}, and
#'  \code{\link{bokeh}}. Otherwise, plots should resemble those
#'  made using \code{perb(rbind(x@counts, y@counts))}.
#'
#' @param x,y A \code{propr} object.
#' @inheritParams perb
#'
#' @return An abstracted \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{prop2prob}}
#'
#' @examples
#' library(propr)
#' data(mail)
#' mail1 <- mail[1:2, ]
#' mail2 <- mail[3:4, ]
#' rho1 <- perb(mail1)
#' rho2 <- perb(mail2)
#' abstract(rho1, rho2)
#' @export
abstract <- function(x, y, select){

  differentialCheck(x, y, forceBoth = TRUE)

  if(!missing(select)){

    x <- subset(x, select = select)
    y <- subset(y, select = select)
  }

  rho <- new("propr")
  rho@counts <- rbind(x@counts, y@counts)
  rho@logratio <- rbind(x@logratio, y@logratio) # clr(x) transforms subject vectors
  rho@matrix <- tanh(atanh(x@matrix) - atanh(y@matrix))
  diag(rho@matrix) <- 1

  return(rho)
}
