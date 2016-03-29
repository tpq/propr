#' @describeIn propr Method to show \code{propr} object.
#'
#' @param object,x An object of class \code{propr}.
#' @export
setMethod("show", "propr",
          function(object){

            cat("@counts summary:",
                nrow(object@counts), "features by", ncol(object@counts), "subjects\n")

            cat("@logratio summary:",
                nrow(object@logratio), "features by", ncol(object@logratio), "subjects\n")

            cat("@matrix summary:",
                nrow(object@matrix), "features by", ncol(object@matrix), "features\n")

            cat("@pairs summary:",
                nrow(object@pairs), "feature pairs\n")
          }
)

#' @describeIn propr Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param subset Subsets via \code{object@counts[subset, ]}.
#' @export
setMethod("subset", signature(x = "propr"),
          function(x, subset){

            x@counts <- x@counts[subset, ]
            x@logratio <- x@logratio[subset, ]
            x@matrix <- x@matrix[subset, subset]
            x@pairs <- proprPairs(x@matrix)

            return(x)
          }
)

#' @describeIn propr Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param i,j,drop Subsets via \code{object@pairs[i, j, drop]}.
#' @export
setMethod('[', signature(x = "propr"),
          function(x, i, j, drop){

            if(!missing(j)){

              return(x@pairs[i, j, drop])

            }else{

              x@pairs <- x@pairs[i, j, drop]
              index <- unique(c(x@pairs$feature1, x@pairs$feature2))
              x@matrix <- x@matrix[index, index]
              x@logratio <- x@logratio[index, ]
              x@counts <- x@counts[index, ]

              return(x)
            }
          }
)

#' @describeIn propr Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param name Subsets via \code{object@pairs[, name]}.
#' @export
setMethod('$', signature(x = "propr"),
          function(x, name){

            return(x@pairs[, name])
          }
)
