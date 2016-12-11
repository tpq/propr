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
#' @inheritParams slate
#'
#' @return Returns a \code{data.table} of p-values.
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

  if(utils::object.size(x) > 10^9 & prompt){

    message("Uh oh! A large proportionality matrix was detected (> 1 GB).\n",
            "This operation requires about 5 times the size of 'x' in additional RAM.\n",
            "Are you sure you want to calculate ",
            ncol(x@matrix) * (ncol(x@matrix)-1) * 1/2, " p-values?\n",
            "0: Nevermind\n1: Proceed\n2: Hmm...")

    response <- readline(prompt = "Which do you choose? ")
    if(!response == 1) stop("prop2prob method aborted.")
  }

  differentialCheck(x, y, forceBoth = FALSE)

  X <- linRcpp(x@matrix, x@logratio[])
  z <- lltRcpp(X)
  var <- urtRcpp(X)
  rm(X); gc()

  if(!missing(y)){

    Y <- linRcpp(y@matrix, y@logratio[])
    z <- z - lltRcpp(Y)
    var <- var + urtRcpp(Y)
    rm(Y); gc()
  }

  # Calculate normal deviate
  z <- suppressWarnings(abs(z / sqrt(var)))
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

  ind <- is.na(dt$Probability)
  if(any(ind)){

    message("Alert: Removing NAs due to NaN variance or alr-transformation.")
    return(dt[!ind, ])

  }else{

    return(dt)
  }
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
#'  \code{tanh(atanh(x@matrix) - atanh(y@matrix))}. The \code{@@pairs}
#'  slot contains an index of all statistically significant pairs,
#'  toggled via the argument \code{dt}.
#'
#' Visualizing the difference matrix with \code{dendrogram} may
#'  help summarize the results of \code{prop2prob}. Note that the
#'  difference matrix now also informs co-cluster assignment for
#'  the \code{bucket}, \code{prism}, and \code{bokeh} plots.
#'  Otherwise, most abstracted plots should match those made using
#'  \code{perb(rbind(x@counts, y@counts))}.
#'
#' @param x,y A \code{propr} object.
#' @param dt A \code{data.table}. The result from
#'  \code{prop2prob(x, y)}.
#' @param colBy A character string. The column in \code{dt}
#'  used to select statistically significant pairs.
#' @param cutoff A numeric scalar. The value of \code{colBy}
#'  used to select statistically significant pairs.
#'
#' @return Returns an abstracted \code{propr} object.
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
abstract <- function(x, y, dt, colBy = "Adjusted", cutoff = .01){

  differentialCheck(x, y, forceBoth = TRUE)

  if(!missing(dt)){

    if(!identical(colnames(dt), c("Partner", "Pair", "Probability", "Adjusted"))){
      stop("Uh oh! Provided 'dt' object does not have the expected column names.")
    }

    message("Alert: Indexing statistically significant pairs.")
    small <- dt[dt[, colBy] < cutoff, ]
    if(nrow(small) == 0) stop("No pairs pass the selected cutoff criteria.")

    x@pairs <- coordToIndex(small$Partner, small$Pair, nrow(x@matrix))
    x <- simplify(x)

    y@pairs <- coordToIndex(small$Partner, small$Pair, nrow(y@matrix))
    y <- simplify(y)
  }

  rho <- new("propr")
  rho@counts <- rbind(x@counts, y@counts)
  rho@logratio <- rbind(x@logratio, y@logratio)
  rho@matrix <- tanh(atanh(x@matrix) - atanh(y@matrix))
  rho@pairs <- union(x@pairs, y@pairs)
  diag(rho@matrix) <- 1

  ind <- is.na(rho@matrix)
  if(any(ind)){

    message("Alert: Replacing NAs due to NaN variance with 0.")
    rho@matrix[ind] <- 0
  }

  return(rho)
}
