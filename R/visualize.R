#' Visualize Proportionality
#'
#' @description
#' \code{smear:}
#'  Plots log-ratio transformed abundances pairwise.
#'  Index-aware, meaning that it only plots pairs indexed
#'  in \code{@@pairs}, unless no pairs are indexed.
#'
#' \code{dendrogram:}
#'  Plots a clustering of proportionality matrix.
#'  Index-aware, meaning that it only plots pairs indexed
#'  in \code{@@pairs}, unless no pairs are indexed.
#'  Heatmap intensity is not scaled.
#'
#' \code{slate:}
#'  Builds a table of VLR, VLS, and proportionality
#'  for each feature pair in a \code{propr} object. If the
#'  argument \code{k} is provided, the table will also
#'  include co-cluster membership.
#'
#' \code{prism:}
#'  Plots the variance of the ratio of the log-ratio transformed
#'  feature pair (VLR) versus the sum of the individual variances
#'  of each log-ratio transformed feature (VLS). The ratio of
#'  the VLR to the VLS equals \code{1 - rho}. As such, we use
#'  here seven rainbow colored lines to indicate where \code{rho}
#'  equals \code{[.01, .05, .50, 0, 1.50, 1.95, 1.99]}, going
#'  from red to violet.
#'
#' \code{bokeh:}
#'  Plots the feature variances for each log-ratio transformed
#'  feature pair in the \code{propr} object. Highly proportional
#'  pairs will aggregate near the \code{y = x} diagonal.
#'  Clusters that appear toward the top-right of the
#'  figure contain features with highly variable abundance across
#'  all samples. Clusters that appear toward the
#'  bottom-left of the figure contain features with fixed
#'  abundance across all samples. Uses a log scale.
#'
#' \code{bucket:}
#'  Plots an estimation of the degree to which a feature pair
#'  differentiates the experimental groups versus the
#'  measure of the proportionality between that feature pair.
#'  "It's pronounced, 'bouquet'." - Hyacinth Bucket
#'
#' \code{mds:}
#'  Plots the first two principal components as calculated
#'  using the log-ratio transformed feature vectors. This
#'  provides a statistically valid alternative to
#'  conventional principal components analysis (PCA) used
#'  in multi-dimensional scaling (MDS) plotting.
#'  For more information, see <DOI:10.1139/cjm-2015-0821>.
#'
#' \code{snapshot:}
#'  Plots the log-ratio transformed feature abundance as
#'  a heatmap, along with the respective dendrograms.
#'  Heatmap intensity is not scaled.
#'
#' \code{cytescape:}
#'  Builds a table of indexed pairs and proportionality.
#'  In doing so, this function displays a preview of the
#'  interaction network, built using \code{igraph}.
#'  We recommend using the result as input to a
#'  network visualization tool like Cytoscape.
#'
#' @return
#' \code{smear, mds:}
#'  Returns a \code{ggplot} object.
#'
#' \code{dendrogram, snapshot:}
#'  Returns a \code{dendrogram} object.
#'
#' \code{slate:}
#'  Returns a \code{data.frame} of all pairwise relationships.
#'  If the argument \code{k} is provided, returns a list of
#'  the \code{data.frame} of pairwise relationships and the
#'  cluster membership.
#'
#' \code{prism, bokeh, bucket:}
#'  Returns cluster membership if \code{k} is provided.
#'  Otherwise, returns a \code{ggplot} object.
#'
#' \code{cytescape:}
#'  Returns a \code{data.frame} of indexed pairs.
#'
#' @param object,rho A \code{propr} object created from \code{perb}.
#'  However, \code{smear}, \code{dendrogram}, and \code{cytescape} will
#'  also accommodate results from \code{phit} and \code{phis}.
#' @param group A character vector. Group or sub-group memberships,
#'  ordered according to the row names in \code{@@counts} and
#'  \code{@@logratio}. Required parameter for \code{bucket}
#'  and optional parameter for \code{mds}.
#' @param k A numeric scalar. The number of clusters. Optional
#'  parameter for \code{bucket}, \code{prism}, and \code{bokeh}.
#'  Providing the argument \code{k} will color feature pairs
#'  by co-cluster membership. In other words, a feature pair
#'  will receive a color if and only if both features belong
#'  to same the cluster (calculated using \code{hclust}).
#' @param prompt A logical scalar. Set to \code{FALSE} to disable
#'  the courtesy prompt when working with big data.
#' @param plotly A logical scalar. Set to \code{TRUE} to produce
#'  a dynamic plot using the \code{plotly} package.
#' @param minPairs An integer scalar. Subsets the interaction
#'  network to exclude any pair without a node that participates
#'  in at least this many total pairs. Required parameter
#'  for \code{cytescape}.
#'
#' @importFrom stats var as.dist as.formula lm aov cutree prcomp dist
#' @name visualize
NULL

#' @rdname visualize
#' @export
smear <- function(rho, prompt = TRUE, plotly = FALSE){

  rho <- plotCheck(rho, prompt = FALSE, plotly = plotly, indexNaive = FALSE)

  if(length(rho@pairs) == 0){

    message("Alert: Generating plot using all feature pairs.")
    V <- indexPairs(rho@matrix, "all")

  }else{

    message("Alert: Generating plot using indexed feature pairs.")
    V <- rho@pairs
  }

  # Build coordinates from index
  coord <- indexToCoord(V, nrow(rho@matrix))
  i.feat <- sort(union(coord[[1]], coord[[2]]))
  if(prompt) promptCheck(length(i.feat))

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

#' @rdname visualize
#' @export
dendrogram <- function(rho, prompt = TRUE, plotly = FALSE){

  dendroCheck()
  rho <- plotCheck(rho, prompt = FALSE, plotly = plotly, indexNaive = FALSE)

  if(length(rho@pairs) == 0){

    message("Alert: Generating plot using all feature pairs.")
    V <- indexPairs(rho@matrix, "all")

  }else{

    message("Alert: Generating plot using indexed feature pairs.")
    V <- rho@pairs
  }

  # Build coordinates from index
  coord <- indexToCoord(V, nrow(rho@matrix))
  i.feat <- sort(union(coord[[1]], coord[[2]]))
  if(prompt) promptCheck(length(i.feat))

  # Make smaller propr object
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

#' @rdname visualize
#' @export
slate <- function(rho, k, prompt = TRUE, plotly = FALSE){

  rho <- plotCheck(rho, prompt = prompt, plotly = plotly, indexNaive = TRUE)

  # Calculate log-ratio transformed variances
  feat.each <- colnames(rho@logratio)
  var.ratio <- vlrRcpp(rho@counts[])
  var.each <- apply(rho@logratio, 2, var)

  # Cluster if k is provided
  if(!missing(k)){

    # Convert rho into dist matrix
    # See reference: http://research.stowers-institute.org/
    #  mcm/efg/R/Visualization/cor-cluster/index.htm
    dist <- as.dist(1 - abs(rho@matrix))
    clust <- cutree(fastcluster::hclust(dist), k = k)
  }

  # Build table components
  llt <- ncol(var.ratio) * (ncol(var.ratio)-1) * 1/2 # size of lower-left triangle
  feat1 <- vector("numeric", llt) # feature name 1
  feat2 <- vector("numeric", llt) # feature name 2
  vl1 <- vector("numeric", llt) # var log feature 1
  vl2 <- vector("numeric", llt) # var log feature 2
  vls <- vector("numeric", llt) # var log sum
  vlr <- vector("numeric", llt) # var log ratio
  rho <- vector("numeric", llt) # 1 - vlr/vls
  col <- vector("numeric", llt) # co-cluster
  count <- 1
  for(j in 2:nrow(var.ratio)){
    for(i in 1:(j-1)){

      feat1[count] <- feat.each[j]
      feat2[count] <- feat.each[i]
      vl1[count] <- var.each[j]
      vl2[count] <- var.each[i]
      vls[count] <- var.each[j] + var.each[i]
      vlr[count] <- var.ratio[j, i]
      rho[count] <- 1 - vlr[count] / vls[count]

      # Since each col initializes as zero
      if(!missing(k)) if(clust[i] == clust[j]) col[count] <- clust[i]

      count <- count + 1
    }
  }

  final <- data.frame("Partner" = feat1, "Pair" = feat2,
                      "VL1" = vl1, "VL2" = vl2,
                      "VLR" = vlr, "VLS" = vls,
                      "rho" = rho)

  if(!missing(k)){

    final$CoCluster <- as.character(col)
    return(list(final, clust))

  }else{

    return(final)
  }
}

#' @rdname visualize
#' @export
bucket <- function(rho, group, k, prompt = TRUE, plotly = FALSE){ # pronounced bouquet

  df <- slate(rho, k, prompt, plotly)

  if(!missing(k)){

    clust <- df[[2]]
    df <- df[[1]]

  }else{

    df$CoCluster <- as.character(0)
  }

  # Calculate discriminating power of each feature
  numfeats <- ncol(rho@logratio)
  data <- data.frame(rho@logratio, group)
  p.val <- vector("numeric", numfeats)
  for(i in 1:numfeats){

    formula <- as.formula(paste0(colnames(data)[i], "~group"))
    fit <- lm(formula, data)
    res <- summary(aov(fit))
    p.val[i] <- res[[1]]$'Pr(>F)'[1]
  }

  # Add discriminating power to slate result
  llt <- ncol(rho@matrix) * (ncol(rho@matrix)-1) * 1/2
  df$Score <- 0
  count <- 1
  for(j in 2:numfeats){
    for(i in 1:(j-1)){
      df$Score[count] <- -log(p.val[i] * p.val[j])
      count <- count + 1
    }
  }

  g <-
    ggplot2::ggplot(
      df, ggplot2::aes_string(x = "rho", y = "Score", Partner = "Partner", Pair = "Pair")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "CoCluster")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::xlab("Proportionality Between Features (rho)") +
    ggplot2::ylab("Discrimination Between Groups") +
    ggplot2::ggtitle("Group Discrimination by Feature Pair") +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(0, max(df$Score)) +
    ggplot2::geom_hline(yintercept = -log(.05 / nrow(df)), color = "lightgrey") +
    ggplot2::geom_hline(yintercept = -log(.05^2 / nrow(df)), color = "black")

  if(plotly){

    return(plotly::ggplotly(g))

  }else{

    plot(g)

    if(!missing(k)){

      return(clust)
    }
  }

  return(g)
}

#' @rdname visualize
#' @export
prism <- function(rho, k, prompt = TRUE, plotly = FALSE){

  df <- slate(rho, k, prompt, plotly)

  if(!missing(k)){

    clust <- df[[2]]
    df <- df[[1]]

  }else{

    df$CoCluster <- as.character(0)
  }

  g <-
    ggplot2::ggplot(df, ggplot2::aes_string(x = "VLS", y = "VLR", rho = "rho",
                                            Partner = "Partner", Pair = "Pair")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "CoCluster")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::xlab("Variance of the Log Sum (VLS)") +
    ggplot2::ylab("Variance of the Log Ratio (VLR)") +
    ggplot2::ggtitle("Distribution of *lr-transformed Variance") +
    ggplot2::xlim(0, max(df$VLS)) +
    ggplot2::ylim(0, max(df$VLR)) +
    ggplot2::geom_abline(slope = 1.99, intercept = 0, color = "violet") +
    ggplot2::geom_abline(slope = 1.95, intercept = 0, color = "purple") +
    ggplot2::geom_abline(slope = 1.50, intercept = 0, color = "blue") +
    ggplot2::geom_abline(slope = 1.00, intercept = 0, color = "green") +
    ggplot2::geom_abline(slope = 0.50, intercept = 0, color = "yellow") +
    ggplot2::geom_abline(slope = 0.05, intercept = 0, color = "orange") +
    ggplot2::geom_abline(slope = 0.01, intercept = 0, color = "red")

  if(plotly){

    return(plotly::ggplotly(g))

  }else{

    plot(g)

    if(!missing(k)){

      return(clust)
    }
  }

  return(g)
}

#' @rdname visualize
#' @export
bokeh <- function(rho, k, prompt = TRUE, plotly = FALSE){

  df <- slate(rho, k, prompt, plotly)

  if(!missing(k)){

    clust <- df[[2]]
    df <- df[[1]]

  }else{

    df$CoCluster <- as.character(0)
  }

  df$logVL1 <- log(df$VL1)
  df$logVL2 <- log(df$VL2)

  g <-
    ggplot2::ggplot(
      df, ggplot2::aes_string(x = "logVL1", y = "logVL2", VLS = "VLS", VLR = "VLR",
                              Partner = "Partner", Pair = "Pair")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "CoCluster", alpha = "rho")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::scale_alpha_continuous(limits = c(-1, 1), name = "Proportionality") +
    ggplot2::xlab("Log *lr-transformed Variance[1]") +
    ggplot2::ylab("Log *lr-transformed Variance[2]") +
    ggplot2::ggtitle("Distribution of *lr-transformed Variance") +
    ggplot2::xlim(min(c(df$logVL1, df$logVL2)), max(c(df$logVL1, df$logVL2))) +
    ggplot2::ylim(min(c(df$logVL1, df$logVL2)), max(c(df$logVL1, df$logVL2))) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "lightgrey")

  if(plotly){

    return(plotly::ggplotly(g))

  }else{

    plot(g)

    if(!missing(k)){

      return(clust)
    }
  }

  return(g)
}

#' @rdname visualize
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

#' @rdname visualize
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

#' @rdname visualize
#' @export
cytescape <- function(object, minPairs = 2){

  packageCheck("igraph")

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  # Prepare data
  rho <- object@matrix[object@pairs]
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  df <- data.frame("Partner" = coords[[1]], "Pair" = coords[[2]], rho)

  # Remove extraneous pairs
  keep <- which(table(c(df$Partner, df$Pair)) >= minPairs)
  sub <- df[df$Partner %in% keep | df$Pair %in% keep, ]

  # Build and color igraph
  if(nrow(sub) == 0) stop("No pairs remain after filter.")
  g <- igraph::graph_from_data_frame(sub, directed = FALSE)
  igraph::V(g)$color <- "white"
  colors <- rep("yellow", nrow(sub))
  colors[sub$rho > .25] <- "orange"
  colors[sub$rho > .75] <- "red"
  colors[sub$rho < -.25] <- "green"
  colors[sub$rho < -.75] <- "blue"
  igraph::E(g)$color <- colors
  plot(g, vertex.size = 2, vertex.label = NA)

  # Retrieve node names
  names <- colnames(object@logratio)
  if(!is.null(names)){
    sub$Partner <- names[sub$Partner]
    sub$Pair <- names[sub$Pair]
  }

  return(sub)
}
