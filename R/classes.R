#' @useDynLib propr
#' @importFrom Rcpp sourceCpp
NULL

#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about how to calculate metrics of proportionality,
#'  see the help file for the functions \code{\link{phit}}
#'  and \code{\link{perb}}.
#'
#' To learn more about the resultant \code{propr} class object, see
#'  the help file for the class definition \code{\link{propr-class}}.
#'
#' To learn more about compositional data analysis, and its relevance
#'  to biological count data, see the bundled vignette.
#'
#' To learn more about \code{propr} class methods, see below.
#'
#' @name propr
NULL

#' An S4 class to hold results from proportionality analysis.
#'
#' @slot counts A matrix. Stores the original "count matrix" input.
#' @slot logratio A matrix. Stores the log-ratio transformed "count matrix".
#' @slot matrix A matrix. Stores the proportionality matrix calculated by
#'  \code{phiRcpp} or \code{rhoRcpp}.
#' @slot pairs A vector. Indexes the proportionality metrics of interest.
#'
#' @seealso \code{\link{propr}}, \code{\link{phit}}, \code{\link{perb}}
#'
#' @examples
#' set.seed(12345)
#' N <- 100
#' X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
#'                 c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
#' rho <- perb(X, ivar = 0)
#' rho[">", .99]
#' subset(rho, 1:2)
#' @export
setClass("propr",
         slots = c(
           counts = "matrix",
           logratio = "matrix",
           matrix = "matrix",
           pairs = "numeric"
         )
)
