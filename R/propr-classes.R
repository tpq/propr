#' @name propr
#' @rdname propr
NULL

#' @rdname propr
#' @section Class:
#'
#' An S4 class to hold results from proportionality analysis.
#'
#' @slot counts A data.frame. Stores the original "counts matrix" input.
#' @slot logratio A data.frame. Stores the log-ratio transformed "counts matrix".
#' @slot matrix A matrix. Stores the proportionality matrix calculated by \code{phit} or \code{perb}.
#' @slot pairs A data.frame. Projects the proportionality matrix pairwise.
#'
#' @seealso \code{\link{phit}}, \code{\link{perb}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 2000 * 22, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 2000, ncol = 22)
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
