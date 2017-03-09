#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about how to calculate proportionality, see
#'  \code{\link{proportionality}}.
#'
#' To learn more about proportionality plots, see
#'  \code{\link{visualize}}.
#'
#' To learn more about differential proportionality, see
#'  \code{\link{prop2prob}}.
#'
#' To learn more about compositional data analysis, and its relevance
#'  to biological count data, see the bundled vignette.
#'
#' @slot counts A matrix. Stores the original "count matrix" input.
#' @slot logratio A matrix. Stores the log-ratio transformed "count matrix".
#' @slot matrix A matrix. Stores the proportionality matrix calculated by
#'  \code{phiRcpp} or \code{rhoRcpp}.
#' @slot pairs A vector. Indexes the proportionality metrics of interest.
#'
#' @param object,x An object of class \code{propr}.
#' @param subset Subsets via \code{object@counts[subset, ]}.
#'  Use this argument to rearrange subject order.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Use this argument to rearrange feature order.
#' @param i Operation used for the subset indexing. Select from
#'  "==", "=", ">", ">=", "<", "<=", "!=", or "all".
#' @param j Reference used for the subset indexing. Provide a numeric
#'  value to which to compare the proportionality measures in the
#'  \code{@@matrix} slot.
#' @param tiny A logical scalar. Toggles whether to pass the indexed
#'  result through \code{\link{simplify}}.
#' @param y Missing. Ignore. Leftover from the generic method
#'  definition.
#' @inheritParams visualize
#'
#' @name propr
#' @useDynLib propr, .registration = TRUE
#' @importFrom methods show new
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname propr
#' @export
setClass("propr",
         slots = c(
           counts = "matrix",
           logratio = "matrix",
           matrix = "matrix",
           pairs = "numeric"
         )
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propr} object.
#' @export
setMethod("show", "propr",
          function(object){

            cat("@counts summary:",
                nrow(object@counts), "subjects by", ncol(object@counts), "features\n")

            cat("@logratio summary:",
                nrow(object@logratio), "subjects by", ncol(object@logratio), "features\n")

            cat("@matrix summary:",
                nrow(object@matrix), "features by", ncol(object@matrix), "features\n")

            if(length(object@pairs) > 0 | nrow(object@matrix) == 0){

              cat("@pairs summary:", length(object@pairs), "feature pairs\n")

            }else{

              cat("@pairs summary: index with `[` method\n")
            }
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{subset:} Method to subset \code{propr} object.
#' @export
setMethod("subset", signature(x = "propr"),
          function(x, subset, select){

            if(missing(subset)) subset <- 1:nrow(x@counts)
            if(missing(select)) select <- 1:ncol(x@counts)

            if(is.character(select)){

              select <- match(select, colnames(x@counts))
            }

            x@counts <- x@counts[subset, select, drop = FALSE]
            x@logratio <- x@logratio[subset, select, drop = FALSE]
            x@matrix <- x@matrix[select, select, drop = FALSE]

            if(length(x@pairs) > 0){

              message("Alert: User must repopulate @pairs slot after `subset`.")
              x@pairs <- vector("numeric")
            }

            return(x)
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{[:} Method to subset \code{propr} object.
#' @aliases [,propr-method
#' @docType methods
#' @export
setMethod('[', signature(x = "propr", i = "ANY", j = "ANY"),
          function(x, i = "all", j, tiny = FALSE){

            if(i == "all"){

              x@pairs <- indexPairs(x@matrix, "all")
              return(x)
            }

            if(!i %in% c("==", "=", ">", ">=", "<", "<=", "!=")){

              stop("Operator not recognized. Index using e.g., `prop[\">\", .95]`.")
            }

            if(missing(j) | !is.numeric(j) | length(j) != 1){

              stop("Reference not found. Index using e.g., `prop[\">\", .95]`.")
            }

            newPairs <- indexPairs(x@matrix, i, j)
            if(length(newPairs) == 0) message("Alert: Method failed to index any pairs.")
            if(length(x@pairs) > 0) message("Alert: Appending prior index.")
            x@pairs <- sort(union(x@pairs, newPairs))

            if(tiny){

              return(simplify(x))

            }else{

              return(x)
            }
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{plot:} Method to plot \code{propr} object.
#' @export
setMethod("plot", signature(x = "propr", y = "missing"),
          function(x, y, prompt = TRUE, plotly = FALSE){

            smear(x, prompt = prompt, plotly = plotly)
          }
)

#' @rdname propr
#' @section Functions:
#' \code{simplify:}
#'  This convenience function takes an indexed \code{propr} object
#'  and subsets the object based on that index. Then, it populates the
#'  \code{@@pairs} slot of the new object with an updated version
#'  of the original index. You can call \code{simplify} from within the
#'  \code{[} method using the argument \code{tiny}.
#' @export
simplify <- function(object){

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  # Subset propr object based on index
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  selection <- sort(union(coords[[1]], coords[[2]]))
  object@pairs <- vector("numeric")
  new <- subset(object, select = selection)

  # Repopulate the pairs slot
  for(i in 1:length(coords[[1]])){

    coords[[1]][i] <- which(selection == coords[[1]][i])
    coords[[2]][i] <- which(selection == coords[[2]][i])
  }

  new@pairs <- (coords[[2]] - 1) * nrow(new@matrix) + (coords[[1]] - 1) + 1

  return(new)
}

#' @rdname propr
#' @section Functions:
#' \code{adjacent:}
#'  This function uses pairs indexed in the \code{@@pairs}
#'  slot to build a symmetric adjacency matrix.
#' @export
adjacent <- function(object){

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  N <- nrow(object@matrix)
  mat <- matrix(0, N, N)
  mat[object@pairs] <- 1
  diag(mat) <- 1
  symRcpp(mat)

  adj <- object
  adj@matrix <- mat

  return(adj)
}
