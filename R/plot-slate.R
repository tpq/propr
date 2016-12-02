#' Build \code{propr} Table
#'
#' This function builds a table of VLR, VLS, and rho
#'  for each feature pair in a \code{propr} object. If the
#'  argument \code{k} is provided, the table will also
#'  include co-cluster membership (as used in the
#'  \code{\link{bucket}}, \code{\link{prism}}, and
#'  \code{\link{bokeh}} plots).
#'
#' @inheritParams bucket
#' @return Returns a \code{data.frame} of pairwise relationships.
#'  If the argument \code{k} is provided, returns a list of
#'  the \code{data.frame} of pairwise relationships and the
#'  cluster membership.
#'
#' @importFrom stats var
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

#' Make Bucket Plot
#'
#' Plots an estimation of the degree to which a feature pair
#'  differentiates the experimental groups versus the
#'  measure of the proportionality between that feature pair.
#'
#' Providing the argument \code{k} will color feature pairs
#'  by co-cluster membership. In other words, a feature pair
#'  will receive a color if and only if both features belong
#'  to same the cluster (calculated using \code{hclust}).
#'
#' "It's pronounced, 'bouquet'." - Hyacinth Bucket
#'
#' @param rho A \code{propr} object created from \code{\link{perb}}.
#'  However, \code{\link{smear}} and \code{\link{dendrogram}} will also
#'  accommodate a \code{propr} object created from \code{\link{phit}}.
#' @param group A character vector. Group or sub-group memberships,
#'  ordered according to the column names in \code{@@counts} and
#'  \code{@@logratio}. Required parameter for \code{\link{bucket}}
#'  and optional parameter for \code{\link{mds}}.
#' @param k A numeric scalar. The number of groups into which to
#'  cluster the subjects. Clusters calculated  using \code{hclust}.
#'  Optional parameter for \code{\link{bucket}}, \code{\link{prism}},
#'  and \code{\link{bokeh}}.
#' @param prompt A logical scalar. Set to \code{FALSE} to disable
#'  the courtesy prompt when working with big data.
#' @param plotly A logical scalar. Set to \code{TRUE} to produce
#'  a dynamic plot using the \code{plotly} package.
#'
#' @return Returns cluster membership if \code{k} is provided.
#'  Otherwise, returns a \code{ggplot} object.
#'
#' @importFrom stats as.formula lm aov cutree
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

#' Make Prism Plot
#'
#' Plots the variance of the ratio of the log-ratio transformed
#'  feature pair (vlr) versus the sum of the individual variances
#'  of each log-ratio transformed feature (vls). The ratio of
#'  the vlr to the vls equals \code{1 - rho}. As such, we use
#'  here seven rainbow colored lines to indicate where \code{rho}
#'  equals \code{[.01, .05, .50, 0, 1.50, 1.95, 1.99]}, going
#'  from red to violet.
#'
#' Providing the argument \code{k} will color feature pairs
#'  by co-cluster membership. In other words, a feature pair
#'  will receive a color if and only if both features belong
#'  to same the cluster (calculated using \code{hclust}).
#'
#' @inheritParams bucket
#' @return Returns cluster membership if \code{k} is provided.
#'  Otherwise, returns a \code{ggplot} object.
#'
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

#' Make Bokeh Plot
#'
#' Plots the feature variances for each log-ratio transformed
#'  feature pair in the \code{propr} object. Highly proportional
#'  pairs will aggregate near the \code{y = x} diagonal.
#'  Clusters that appear toward the top-right of the
#'  figure contain features with highly variable abundance across
#'  all samples. Clusters that appear toward the
#'  bottom-left of the figure contain features with fixed
#'  abundance across all samples. Uses a log scale.
#'
#' Providing the argument \code{k} will color feature pairs
#'  by co-cluster membership. In other words, a feature pair
#'  will receive a color if and only if both features belong
#'  to same the cluster (calculated using \code{hclust}).
#'
#' @inheritParams bucket
#' @return Returns cluster membership if \code{k} is provided.
#'  Otherwise, returns a \code{ggplot} object.
#'
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
