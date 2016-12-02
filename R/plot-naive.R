#' Make MDS Plot
#'
#' Plots the first two principal components as calculated
#'  using the log-ratio transformed feature vectors. This
#'  provides a statistically valid alternative to
#'  conventional principal components analysis (PCA) used
#'  in multi-dimensional scaling (MDS) plotting.
#'
#' For more information, see DOI:10.1139/cjm-2015-0821.
#'
#' @inheritParams bucket
#' @return Returns a \code{ggplot} object.
#'
#' @importFrom stats prcomp
#' @export
mds <- function(rho, group, prompt = TRUE, plotly = FALSE){

  rho <- plotCheck(rho, prompt = prompt, plotly = plotly, indexNaive = TRUE)

  if(missing(group)){

    group <- rep("None", nrow(rho@logratio))
  }

  df <- data.frame("ID" = rownames(rho@logratio), "Group" = as.character(group),
                   prcomp(rho@logratio)$x[, c(1, 2)])
  g <-
    ggplot2::ggplot(ggplot2::aes_string(ID = "ID"), data = df) +
    ggplot2::geom_point(
      ggplot2::aes_string(x = "PC1", y = "PC2", colour = "Group")) +
    ggplot2::theme_bw() +
    ggplot2::xlab("First Component") + ggplot2::ylab("Second Component") +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Group") +
    ggplot2::ggtitle("*lr-transformed MDS Plot")

  if(plotly){

    return(plotly::ggplotly(g))

  }else{

    g <- g + ggplot2::geom_text(
      ggplot2::aes_string(x = "PC1", y = "PC2", label = "ID", colour = "Group"),
      data = df, size = 3, vjust = -1)
    plot(g)
  }

  return(g)
}

#' Make Snapshot Plot
#'
#' Plots the log-ratio transformed feature abundance as
#'  a heatmap, along with the respective dendrograms.
#'  Heatmap intensity is not scaled.
#'
#' @inheritParams bucket
#' @return Returns a dendrogram object made from \code{hclust}.
#'
#' @importFrom stats dist
#' @export
snapshot <- function(rho, prompt = TRUE, plotly = FALSE){

  dendroCheck()
  rho <- plotCheck(rho, prompt = prompt, plotly = plotly, indexNaive = TRUE)

  # Build a blank figure
  p_empty <- ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::theme_minimal()

  # Build the column and row dendrograms
  dd.col <- stats::as.dendrogram(fastcluster::hclust(dist(rho@logratio)))
  dd.row <- stats::as.dendrogram(fastcluster::hclust(dist(t(rho@logratio))))
  px <- ggdend(dd.row)
  py <- ggdend(dd.col) + ggplot2::coord_flip()

  # Build the heatmap
  col.ord <- stats::order.dendrogram(dd.col)
  row.ord <- stats::order.dendrogram(dd.row)
  xx <- rho@logratio[col.ord, row.ord]
  df <- as.data.frame(xx)
  df$row <- rownames(df)
  df$row <- with(df, factor(row, levels = row, ordered = TRUE))
  mdf <- reshape2::melt(df, id.vars = "row")
  colnames(mdf) <- c("Sample", "Feature", "lrAbundance")
  p <-
    ggplot2::ggplot(mdf, ggplot2::aes_string(x = "Feature", y = "Sample")) +
    ggplot2::xlab("Features") + ggplot2::ylab("Samples") +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "lrAbundance")) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(axis.text = ggplot2::element_blank())

  if(plotly){

    p <- p + ggplot2::scale_fill_distiller("*lr", palette = "Spectral")
    return(plotly::subplot(px, p_empty, p, py, nrows = 2, margin = 0.01))

  }else{

    p <- p + ggplot2::scale_fill_distiller(palette = "Spectral", guide = FALSE)
    multiplot(px, p, p_empty, py, cols = 2)
    return(dd.col)
  }
}
