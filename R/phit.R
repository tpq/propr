#' Calculate proportionality metric phi.
#'
#' \code{phit} returns a \code{propr} object containing a measure of proportionality.
#'
#' Let d represent any number of features measured across n biological replicates
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with n rows and d columns into a proportionality matrix of d rows and d
#' 	columns containing phi measurements for each feature pair. One can think of the resultant
#' 	matrix as equivalent to a distance matrix, except that it has no symmetry unless forced.
#'
#' @param counts A data.frame or matrix. A "count matrix" with subjects as rows and features as columns.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry by reflecting the "lower left triangle".
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{perb}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' phi <- phit(counts, symmetrize = TRUE)
#' @importFrom methods new
#' @export
phit <- function(counts, symmetrize = TRUE){

    cat("Calculating phi from \"count matrix\".\n")
    prop <- new("propr")
    prop@counts <- as.matrix(counts)
    prop@logratio <- clrRcpp(prop@counts[]) # [] forces copy
    prop@matrix <- phiRcpp(prop@counts[], symmetrize) # [] forces copy
    prop@pairs <- vector("numeric")

    return(prop)
}
