#' @rdname propr
#' @section Methods (by generic):
#' \code{plot:} Method to plot \code{propr} object.
#'
#' @inheritParams bucket
#' @param y Missing. Ignore. Leftover from the generic method definition.
#' @export
setMethod("plot", signature(x = "propr", y = "missing"),
          function(x, y, plotly = FALSE){

            smear(x, plotly = plotly)
          }
)

#' Make Smear Plot
#'
#' Plots *lr-transformed abundances for all indexed pairs.
#'
#' Note that only \code{\link{smear}} and \code{\link{dendrogram}}
#'  produce "index-aware" plots. These functions also accommodate
#'  results from either \code{\link{phit}} or \code{\link{perb}}.
#'
#' @inheritParams bucket
#' @return Returns a \code{ggplot} object.
#'
#' @export
smear <- function(rho, plotly = FALSE){

  rho <- plotCheck(rho, prompt = FALSE, plotly = plotly, indexNaive = FALSE)

  if(length(rho@pairs) == 0){

    message("Alert: Generating plot using all feature pairs.")
    V <- indexPairs(rho@matrix, "all")
    coord <- indexToCoord(V, nrow(rho@matrix))

  }else{

    message("Alert: Generating plot using indexed feature pairs.")
    V <- rho@pairs
    coord <- indexToCoord(V, nrow(rho@matrix))
  }

  # Melt *lr counts by feature pairs
  nsubj <- nrow(rho@logratio)
  L <- length(V) * nsubj
  partner <- vector("character", L)
  pair <- vector("character", L)
  feat1 <- vector("numeric", L)
  feat2 <- vector("numeric", L)
  group <- vector("numeric", length(V) * nsubj)
  for(i in 1:length(V)){

    cat("Shaping pair ", i, "...", sep = "")
    i.order <- order(rho@logratio[, coord$feat1[i]])
    partner[((i-1)*nsubj + 1):((i-1)*nsubj + nsubj)] <- colnames(rho@logratio)[coord$feat1[i]]
    pair[((i-1)*nsubj + 1):((i-1)*nsubj + nsubj)] <- colnames(rho@logratio)[coord$feat2[i]]
    feat1[((i-1)*nsubj + 1):((i-1)*nsubj + nsubj)] <- rho@logratio[, coord$feat1[i]][i.order]
    feat2[((i-1)*nsubj + 1):((i-1)*nsubj + nsubj)] <- rho@logratio[, coord$feat2[i]][i.order]
    group[((i-1)*nsubj + 1):((i-1)*nsubj + nsubj)] <- i
  }

  # Plot *lr-Y by *lr-X
  cat("\n")
  df <- data.frame("X" = feat1, "Y" = feat2, "Group" = group, "Partner" = partner, "Pair" = pair)
  df$Group <- factor(df$Group)
  g <-
    ggplot2::ggplot(
      ggplot2::aes_string(x = "X", y = "Y", Partner = "Partner", Pair = "Pair"), data = df) +
    ggplot2::geom_path(ggplot2::aes_string(colour = "Group")) +
    ggplot2::labs(x = "*lr-transformed Abundance[1]",
                  y = "*lr-transformed Abundance[2]") +
    ggplot2::coord_equal(ratio = 1) + ggplot2::theme_bw() +
    ggplot2::ggtitle("Distribution of *lr-transformed Abundance") +
    ggplot2::theme(legend.position = "none")

  if(plotly){

    return(plotly::ggplotly(g))

  }else{

    plot(g)
  }

  return(g)
}

#' Make Dendrogram Plot
#'
#' Plots a dendrogram and heatmap for all indexed pairs.
#'
#' Note that only \code{\link{smear}} and \code{\link{dendrogram}}
#'  produce "index-aware" plots. These functions also accommodate
#'  results from either \code{\link{phit}} or \code{\link{perb}}.
#'
#' @inheritParams bucket
#' @return Returns a dendrogram object made from \code{hclust}.
#'
#' @importFrom stats as.dist
#' @export
dendrogram <- function(rho, plotly = FALSE){

  dendroCheck()
  rho <- plotCheck(rho, prompt = FALSE, plotly = plotly, indexNaive = FALSE)

  if(length(rho@pairs) == 0){

    message("Alert: Generating plot using all feature pairs.")
    i.feat <- 1:nrow(rho@matrix)

  }else{

    message("Alert: Generating plot using indexed feature pairs.")
    V <- rho@pairs
    coord <- indexToCoord(V, nrow(rho@matrix))
    i.feat <- sort(union(coord[[1]], coord[[2]]))
  }

  rho <- suppressMessages(subset(rho, select = i.feat))

  if(rho@matrix[1, 1] == 0){

    # Convert phi into dis matrix
    dis <- as.dist(rho@matrix)

  }else if(rho@matrix[1, 1] == 1){

    # Convert rho into dis matrix
    # See reference: http://research.stowers-institute.org/
    #  mcm/efg/R/Visualization/cor-cluster/index.htm
    dis <- as.dist(1 - abs(rho@matrix))

  }else{

    stop("Matrix not recognized.")
  }

  # Build a blank figure
  p_empty <- ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::theme_minimal()

  # Build the column and row dendrograms
  dd.col <- stats::as.dendrogram(fastcluster::hclust(dis))
  px <- ggdend(dd.col)
  py <- px + ggplot2::coord_flip()

  # Build the heatmap
  col.ord <- stats::order.dendrogram(dd.col)
  xx <- rho@matrix[col.ord, col.ord]
  df <- as.data.frame(xx)
  colnames(df) <- colnames(rho@logratio)[col.ord]
  rownames(df) <- colnames(df)
  df$row <- rownames(df)
  df$row <- with(df, factor(row, levels = row, ordered = TRUE))
  mdf <- reshape2::melt(df, id.vars = "row")
  colnames(mdf) <- c("Sample", "Feature", "rho")
  p <-
    ggplot2::ggplot(mdf, ggplot2::aes_string(x = "Feature", y = "Sample")) +
    ggplot2::xlab("Features") + ggplot2::ylab("Features") +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "rho")) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(axis.text = ggplot2::element_blank())

  if(plotly){

    p <- p + ggplot2::scale_fill_distiller(limits = c(-1, 1), name = "Proportionality",
                                           palette = "Spectral")
    return(plotly::subplot(px, p_empty, p, py, nrows = 2, margin = 0.01))

  }else{

    p <- p + ggplot2::scale_fill_distiller(limits = c(-1, 1), name = "Proportionality",
                                           palette = "Spectral", guide = FALSE)
    multiplot(px, p, p_empty, py, cols = 2)
    return(dd.col)
  }
}
