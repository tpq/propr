#' Create permuted data
#'
#' This function creates p permuted data matrices
#'
#' This function wraps \code{updatePermutes.propr} and
#'  \code{updatePermutes.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param p The number of permutations to perform.
#' @return A \code{propr} or \code{propd} object with the permutes slot updated.
#' @export
updatePermutes <-
  function(object, p, fixseed = FALSE) {
    if (inherits(object, "propr")) {
      updatePermutes.propr(object, p, fixseed)

    } else if (inherits(object, "propd")) {
      updatePermutes.propd(object, p, fixseed)

    } else{
      stop("Provided 'object' not recognized.")
    }
  }

updatePermutes.propr <-
  function(object, p = 100, fixseed = FALSE) {
    # Shuffle row-wise so that per-sample CLR does not change
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- vector("list", p)
    for (ins in 1:p){
      if (fixseed) set.seed(ins)
      permutes[[ins]] <- t(apply(ct, 1, sample))
    }
    object@permutes <- permutes
    return(object)
  }

updatePermutes.propd <-
  function(object, p = 100, fixseed = FALSE) {
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
    for (col in 1:p){
      if (fixseed) set.seed(col)
      permutes[, col] <- sample(1:nrow(ct))
    }
    object@permutes <- permutes
    return(object)
  }
