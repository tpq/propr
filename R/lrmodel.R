#' An S4 class to hold results from proportionality analysis.
#'
#' @slot MARGIN Specifies the \code{MARGIN} that contains the sample
#'  data. If \code{MARGIN = 1}, rows contain sample data.
#'  If \code{MARGIN = 2}, columns contain sample data.
#' @slot replace Stores the smallest non-zero value of the raw data,
#'  used for replacing zeroes in the independent dataset.
#' @slot type Specifies the type of log-ratio transformation used.
#' @slot rule Specifies the denominators used during the
#'  log-ratio transformation. For a clr-transformation, this contains
#'  the geometric means of the samples. For an alr-transformation,
#'  this contains the feature with fixed abundance.
#' @slot data The raw data as transformed by the \code{rule}.
#' @slot raw The raw data used to build the \code{rule}.
#'
#' @seealso \code{\link{modelCLR}}
#'
#' @examples
#' library(propr)
#' data(mail)
#' model <- modelCLR(mail)
#' lr <- predict(model, mail)
#' @export
setClass("lrmodel",
         slots = c(
           MARGIN = "numeric",
           replace = "numeric",
           type = "character",
           rule = "numeric",
           data = "matrix",
           raw = "matrix"
         )
)

#' Encapsulate a clr-Transformation
#'
#' This function encapsulates the rule for a clr-transformation to apply
#'  to an independent dataset (via \code{predict}).
#'
#' @param x A non-negative count matrix.
#' @param MARGIN If \code{MARGIN = 1}, rows contain sample data.
#'  If \code{MARGIN = 2}, columns contain sample data.
#'
#' @seealso \code{\link{lrmodel-class}}
#'
#' @export
modelCLR <- function(x, MARGIN = 1){

  if(class(x) != "matrix") stop("You must provide a matrix.")
  if(any(0 == x)){
    message("Alert: Replacing 0s in \"count matrix\" with next smallest value.")
    x[x == 0] <- unique(sort(x))[2]
  }

  logX <- log(x)
  means <- apply(logX, MARGIN, mean)
  data <- sweep(logX, MARGIN, means, "-")

  new("lrmodel",
      MARGIN = MARGIN,
      replace = min(x),
      type = "clr",
      rule = means,
      data = data,
      raw = x)
}

#' @describeIn lrmodel Log-ratio transform new data.
#'
#' @param object An \code{lrmodel} object.
#' @param newdata A non-negative count matrix.
#'
#' @export
setMethod("predict", "lrmodel",
          function(object, newdata){

            if(class(newdata) != "matrix") stop("You must provide the new data as a matrix.")
            if(any(0 == newdata)){
              message("Alert: Replacing 0s in new data with integer replacement.")
              newdata[newdata == 0] <- object@replace
            }

            sweep(log(newdata), object@MARGIN, object@rule, "-")
          }
)
