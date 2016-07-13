#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about how to calculate metrics of proportionality,
#'  see the help file for the method definitions \code{\link{phit}}
#'  and \code{\link{perb}}.
#'
#' To learn more about the resultant \code{propr} class object, see
#'  the help file for the class definition \code{\link{propr-class}}.
#'
#' To learn more about compositional data analysis, and its relevance
#'  to biological count data, see the attached vignette.
#'
#' To learn more about \code{propr} class methods, see below.
#'
#' @name propr
NULL

#' An S4 class to hold results from proportionality analysis.
#'
#' @slot counts A data.frame. Stores the original "count matrix" input.
#' @slot logratio A data.frame. Stores the log-ratio transformed "count matrix".
#' @slot matrix A matrix. Stores the proportionality matrix calculated by \code{phit} or \code{perb}.
#' @slot pairs A data.frame. Projects the proportionality matrix pairwise.
#'
#' @seealso \code{\link{propr}}, \code{\link{phit}}, \code{\link{perb}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' prop <- perb(counts, ivar = 0, iter = 0)
#' prop[1:5, ]
#' prop$prop
#' prop[1:5, "prop"]
#' subset(prop, 1:5)
#' @export
setClass("propr",
         slots = c(
           counts = "data.frame",
           logratio = "data.frame",
           matrix = "matrix",
           pairs = "data.frame"
         )
)
