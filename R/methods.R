#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about how to calculate proportionality, see
#'  see \code{\link{proportionality}}.
#'
#' To learn more about \code{propr} plots, see \code{\link{smear}},
#'  \code{\link{dendrogram}}, \code{\link{bucket}}, \code{\link{prism}},
#'  \code{\link{bokeh}}, \code{\link{mds}}, and \code{\link{snapshot}}.
#'
#' To learn more about differential proportionality, see
#'  \code{\link{prop2prob}} and \code{\link{abstract}}.
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
#' @inheritParams bucket
#'
#' @name propr
#' @useDynLib propr
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
