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

#' @describeIn propr Method to plot \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param title A character string. A title for the \code{propr} plot.
#' @export
setMethod("plot", signature(x = "propr", y = "missing"),
          function(x, title = "Pairwise Proportionality"){

            if(!requireNamespace("ggplot2", quietly = TRUE)){
              stop("Uh oh! This plot method depends on ggplot2! ",
                   "Try running: install.packages('ggplot2')")
            }

            if(!requireNamespace("ggthemes", quietly = TRUE)){
              stop("Uh oh! This plot method depends on ggthemes! ",
                   "Try running: install.packages('ggthemes')")
            }

            # Melt *lr counts by feature pairs
            pairs <- vector("list", nrow(x@pairs))
            for(i in 1:nrow(x@pairs)){

              cat("Shaping pair", i, "...")
              pairs[[i]] <- data.frame("x.val" = unlist(x@logratio[x$feature1[i], ]),
                                       "y.val" = unlist(x@logratio[x$feature2[i], ]),
                                       "x.id" = x$feature1[i],
                                       "y.id" = x$feature2[i],
                                       "group" = i)

              pairs[[i]] <- pairs[[i]][order(pairs[[i]]$x.val), ]
            }

            # Plot *lr-Y by *lr-X
            df <- do.call(rbind, pairs)
            p <- ggplot2::ggplot(data = df, aes(x = x.val,
                                                y = y.val,
                                                group = factor(group))) +
              ggplot2::geom_path(aes(colour = factor(group))) +
              ggplot2::labs(x = "Exprssion *LR mRNA[1]",
                            y = "Expression *LR mRNA[2]") +
              ggplot2::coord_equal(ratio = 1) +
              ggthemes::theme_base() +
              ggplot2::theme(legend.position = "none") +
              ggplot2::ggtitle(title)
            plot(p)

            return(p)
          }
)
