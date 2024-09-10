#' Create permuted data
#'
#' This function creates p permuted data matrices
#'
#' This function wraps \code{updatePermutes.propr} and
#'  \code{updatePermutes.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param p The number of permutations to perform. Default is 100.
#' @return A \code{propr} or \code{propd} object with the permutes slot updated.
#' @export
updatePermutes <-
  function(object, p=100) {
    if (inherits(object, "propr")) {
      updatePermutes.propr(object, p)

    } else if (inherits(object, "propd")) {
      updatePermutes.propd(object, p)

    } else{
      stop("Provided 'object' not recognized.")
    }
  }

updatePermutes.propr <-
  function(object, p) {
    # Shuffle row-wise so that per-sample CLR does not change
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- vector("list", p)
    for (ins in 1:p){
      permutes[[ins]] <- t(apply(ct, 1, sample))
    }
    object@permutes <- permutes
    return(object)
  }

updatePermutes.propd <-
  function(object, p) {
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
    for (col in 1:p){
      permutes[, col] <- sample(1:nrow(ct))
    }
    object@permutes <- permutes
    return(object)
  }
