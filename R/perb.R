#' Calculate proportionality metric rho.
#'
#' \code{perb} returns a \code{propr} object containing measures of proportionality.
#'
#' Calculates proportionality metric rho described in Lovell 2015 and expounded
#'  in Erb 2016. Uses centered log-ratio transformation of data by default,
#'  but will use additive log-ratio transformation of data if non-zero
#'  \code{ivar} provided.
#'
#' Let d represent any number of features measured across multiple biological replicates n
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with n rows and d columns into a proportionality matrix of d rows and d
#' 	columns containing rho measurements for each feature pair. One can think of the resultant
#' 	matrix as equivalent to a correlation matrix.
#'
#' @param ivar A numeric scalar. Specificies reference feature for additive log-ratio transformation.
#' @inheritParams phit
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{phit}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' prop <- perb(counts, ivar = 0)
#' @importFrom methods new
#' @importFrom stats ecdf p.adjust
#' @export
perb <- function(counts, ivar = 0){

  cat("Calculating rho from \"count matrix\".\n")
  prop <- new("propr")
  prop@counts <- as.matrix(counts)
  if(ivar != 0){ prop@logratio <- alrRcpp(prop@counts, ivar)
  }else{ prop@logratio <- clrRcpp(prop@counts)}
  prop@matrix <- rhoRcpp(prop@counts, ivar)
  prop@pairs <- vector("numeric")

  return(prop)
}
