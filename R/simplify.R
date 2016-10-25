#' Simplify Indexed Matrix
#'
#' This convenience function takes an indexed \code{\link{propr}} object
#'  and subsets the object based on that index. Then, it populates the
#'  \code{@@pairs} slot of the new object with an updated version
#'  of the original index.
#'
#' @inheritParams propr
#' @return Returns a \code{propr} object.
#'
#' @export
simplify <- function(object){

  if(!class(object) == "propr"){

    stop("Uh oh! You can only simplify an indexed 'propr' object.")
  }

  if(length(object@pairs) == 0){

    stop("Uh oh! You can only simplify an indexed 'propr' object.")
  }

  # Call indexToCoord on indexed propr object
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  selection <- sort(union(coords[[1]], coords[[2]]))
  object@pairs <- vector("numeric")

  # Subset propr object based on index
  new <- subset(object, select = selection)

  # Repopulate the pairs slot
  for(i in 1:length(coords[[1]])){

    coords[[1]][i] <- which(selection == coords[[1]][i])
    coords[[2]][i] <- which(selection == coords[[2]][i])
  }

  new@pairs <- (coords[[1]] - 1) * nrow(new@matrix) + (coords[[2]] - 1) + 1

  return(new)
}
