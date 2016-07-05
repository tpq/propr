#' @rdname propr
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propr} object.
#'
#' @param object,x An object of class \code{propr}.
#' @importFrom methods show
#' @export
setMethod("show", "propr",
          function(object){

            cat("@counts summary:",
                nrow(object@counts), "subjects by", ncol(object@counts), "features\n")

            cat("@logratio summary:",
                nrow(object@logratio), "subjects by", ncol(object@logratio), "features\n")

            cat("@matrix summary:",
                nrow(object@matrix), "features by", ncol(object@matrix), "features\n")

            cat("@pairs summary:",
                nrow(object@pairs), "feature pairs\n")
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{subset:} Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param subset Subsets via \code{object@counts[subset, ]}.
#'  Use this argument to rearrange subject order.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Use this argument to rearrange feature order.
#' @export
setMethod("subset", signature(x = "propr"),
          function(x, subset, select){

            if(missing(subset)) subset <- rownames(x@counts)
            if(missing(select)) select <- colnames(x@counts)

            x@counts <- x@counts[subset, select, drop = FALSE]
            x@logratio <- x@logratio[subset, colnames(x@counts), drop = FALSE]
            x@matrix <- x@matrix[colnames(x@counts), colnames(x@counts), drop = FALSE]
            x@pairs <- proprPairs(x@matrix)

            return(x)
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{[:} Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param i,j,drop Subsets via \code{object@pairs[i, j, drop]}.
#' @aliases [,propr-method
#' @docType methods
#' @export
setMethod('[', signature(x = "propr", i = "ANY", j = "ANY"),
          function(x, i, j, drop){

            if(!missing(j)){

              return(x@pairs[i, j, drop = drop])

            }else{

              x@pairs <- x@pairs[i, j, drop = FALSE]
              index <- unique(c(x@pairs$feature1, x@pairs$feature2))
              x@matrix <- x@matrix[index, index, drop = FALSE]
              x@logratio <- x@logratio[, colnames(x@matrix), drop = FALSE]
              x@counts <- x@counts[, colnames(x@matrix), drop = FALSE]

              return(x)
            }
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{$:} Method to subset \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param name Subsets via \code{object@pairs[, name]}.
#' @export
setMethod('$', signature(x = "propr"),
          function(x, name){

            return(x@pairs[, name])
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{plot:} Method to plot \code{propr} object.
#'
# #' @param x An object of class \code{propr}.
#' @param y Missing. Ignore. Leftover from the generic method definition.
#' @param title A character string. A title for the \code{propr} plot.
#' @export
setMethod("plot", signature(x = "propr", y = "missing"),
          function(x, y, title = "Pairwise Proportionality"){

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
              pairs[[i]] <- data.frame("x.val" = unlist(x@logratio[, x$feature1[i]]),
                                       "y.val" = unlist(x@logratio[, x$feature2[i]]),
                                       "x.id" = x$feature1[i],
                                       "y.id" = x$feature2[i],
                                       "group" = i)

              pairs[[i]] <- pairs[[i]][order(pairs[[i]]$x.val), ]
            }

            # Plot *lr-Y by *lr-X
            df <- do.call(rbind, pairs)
            p <- ggplot2::ggplot(data = df,
                                 ggplot2::aes_string(x = "x.val",
                                                     y = "y.val",
                                                     group = "group")) +
              ggplot2::geom_path(ggplot2::aes(colour = factor(df$group))) +
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

#' @rdname propr
#' @section Methods (by generic):
#' \code{image:} Method to plot \code{propr} object.
#'
#' @param cexRow Numeric. Size of x-axis label.
#' @param cexCol Numeric. Size of y-axis label.
# #' @param object An object of class \code{propr}.
# #' @param title A character string. A title for the \code{propr} plot.
#' @export
setMethod("image", signature(x = "propr"),
          function(x, cexRow = 10, cexCol = 10, title = "*LR Transformed Image"){

            if(!requireNamespace("ggplot2", quietly = TRUE)){
              stop("Uh oh! This plot method depends on ggplot2! ",
                   "Try running: install.packages('ggplot2')")
            }

            # Melt *lr counts by feature
            features <- vector("list", ncol(x@logratio))
            for(i in 1:ncol(x@logratio)){

              features[[i]] <- data.frame("feature" = colnames(x@logratio)[i],
                                          "subject" = rownames(x@logratio),
                                          "value" = unlist(x@logratio[, i]),
                                          stringsAsFactors = FALSE)
            }

            # Image *lr counts
            df <- do.call(rbind, features)
            df$feature <- factor(df$feature, levels = unique(df$feature))
            df$subject <- factor(df$subject, levels = unique(df$subject))
            valMin <- floor(min(df$value))
            valMax <- ceiling(max(df$value))

            # Plot *lr for each subject
            p <- ggplot2::ggplot(ggplot2::aes_string(x = "subject",
                                                     y = "feature"), data = df) +
              ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) +
              ggplot2::scale_fill_gradient2(name = "Scaled *LR Expression",
                                            # low = "grey0", high = "grey70", mid = "grey35",
                                            low = "yellow", high = "red", mid = "orange",
                                            midpoint = 0,
                                            limit = c(valMin, valMax),
                                            breaks = seq(valMin, valMax)) +
              ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             axis.title.y = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_text(size = cexRow),
                             axis.text.y = ggplot2::element_text(size = cexCol),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             legend.position = "bottom") +
              ggplot2::ggtitle(title)
            plot(p)

            return(p)
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{plot:} Method to plot \code{propr} object.
#'
# #' @param object An object of class \code{propr}.
# #' @param title A character string. A title for the \code{propr} plot.
#' @param group A character or numeric vector. Supply feature groups for coloring.
#'  Feature groups expected in the order they appear in \code{@@counts}.
#' @importFrom stats as.dist as.dendrogram hclust order.dendrogram
#' @importFrom grDevices rainbow
#' @export
dendrogram <- function(object, title = "Proportional Clusters", group){

  if(!requireNamespace("dendextend", quietly = TRUE)){
    stop("Uh oh! This plot method depends on dendextend! ",
         "Try running: install.packages('dendextend')")
  }

  if(object@matrix[1, 1] == 0){

    # Convert phi into dist matrix
    dist <- as.dist(object@matrix)

  }else if(object@matrix[1, 1] == 1){

    # Convert rho into dist matrix
    # See reference: http://research.stowers-institute.org/
    #  mcm/efg/R/Visualization/cor-cluster/index.htm
    dist <- as.dist(1 - abs(object@matrix))
  }

  # Align features with groups in data.frame
  if(missing(group)) group <- 1
  colorKey <- data.frame("feature" = colnames(object@matrix),
                         "group" = group, "color" = NA,
                         stringsAsFactors = FALSE)

  # Assign 'n' colors based on 'n' groups
  grps <- unique(colorKey$group)
  colors <- rainbow(length(grps))
  for(i in 1:length(grps)){

    colorKey[colorKey$group == grps[i], "color"] <- colors[i]
  }

  # Build tree and color branches
  dend <- as.dendrogram(hclust(dist))
  dendextend::labels_colors(dend) <- colorKey$color[order.dendrogram(dend)]
  plot(dend, main = title)

  return(dend)
}
