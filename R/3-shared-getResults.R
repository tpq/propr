#' Get Results from Object
#'
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propr} or \code{propd} object.
#'
#' @param object A \code{propr} or \code{propd} object.
#'
#' @return A \code{data.frame} of results.
#'
#' @export
getResults <-
  function(object) {
    df <- object@results
    names <- colnames(object@counts)
    df$Partner <- names[df$Partner]
    df$Pair <- names[df$Pair]
    return(df)
  }
