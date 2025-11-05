#' Create permuted data
#'
#' This function creates p permuted data matrices
#'
#' This function wraps \code{updatePermutes.propr} and
#'  \code{updatePermutes.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param p The number of permutations to perform. Default is 100.
#' @param permutation_option A character string indicating if permute the data
#'  sample-wise or feature-wise. Default is "feature-wise". Note that this flag
#'  is only relevant for \code{propr} objects.
#' @return A \code{propr} or \code{propd} object with the permutes slot updated.
#' @export
updatePermutes <- function(object, p=100, permutation_option = c("feature-wise", "sample-wise")) {
    if (inherits(object, "propr")) {
      updatePermutes.propr(object, p, permutation_option)

    } else if (inherits(object, "propd")) {
      updatePermutes.propd(object, p)

    } else{
      stop("Provided 'object' not recognized.")
    }
  }

updatePermutes.propr <- function(object, p, permutation_option = c("feature-wise", "sample-wise")) {
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- vector("list", p)
    for (ins in 1:p){
      if (permutation_option[1] == "feature-wise") {
        # Permute features
        permutes[[ins]] <- apply(ct, 2, sample)
      } else if (permutation_option[1] == "sample-wise") {
        # Permute samples
        permutes[[ins]] <- t(apply(ct, 1, sample))
      } else {
        stop("Invalid permutation option. Choose either 'feature-wise' or 'sample-wise'.")
      }
    }
    object@permutes <- permutes
    return(object)
  }

updatePermutes.propd <- function(object, p) {
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
    for (col in 1:p){
      permutes[, col] <- sample(1:nrow(ct))
    }
    object@permutes <- permutes
    return(object)
  }