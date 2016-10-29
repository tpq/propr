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
#' @param group A character vector. Group or sub-group memberships,
#'  ordered according to the column names in \code{@@counts} and
#'  \code{@@logratio}. Parameter required for \code{\link{bucket}}
#'  and optional for \code{\link{mds}}.
#' @param k A numeric scalar. The number of groups into which to
#'  cluster the subjects. Clusters calculated  using \code{hclust}.
#'  Optional parameter for \code{\link{bucket}}, \code{\link{prism}},
#'  and \code{\link{bokeh}}.
#' @param prompt A logical scalar. Set to \code{FALSE} to disable
#'  the prompt when plotting a large number of features.
#'
#' @return Returns cluster membership if \code{k} is provided.
#'
#' @importFrom stats as.formula lm aov cutree
#' @export
bucket <- function(rho, group, k, prompt = TRUE){ # pronounced bouquet

  if(suppressWarnings(!requireNamespace("fastcluster", quietly = TRUE))){
    stop("Uh oh! This plot method depends on fastcluster! ",
         "Try running: install.packages('fastcluster')")
  }

  plotCheck(rho, prompt)

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

  # Cluster if k is provided
  if(!missing(k)){

    # Convert rho into dist matrix
    # See reference: http://research.stowers-institute.org/
    #  mcm/efg/R/Visualization/cor-cluster/index.htm
    dist <- as.dist(1 - abs(rho@matrix))
    clust <- cutree(fastcluster::hclust(dist), k = k)
  }

  # Build graph components
  llt <- ncol(rho@matrix) * (ncol(rho@matrix)-1) * 1/2 # size of lower-left triangle
  prop <- vector("numeric", llt)
  weight <- vector("numeric", llt)
  col <- vector("numeric", llt)
  count <- 1
  for(j in 2:numfeats){
    for(i in 1:(j-1)){

      prop[count] <- rho@matrix[j, i]
      weight[count] <- -log(p.val[i] * p.val[j])

      # Since each col initializes as zero
      if(!missing(k)) if(clust[i] == clust[j]) col[count] <- clust[i]

      count <- count + 1
    }
  }

  df <- data.frame(prop, weight, "col" = as.character(col))
  g <-
    ggplot2::ggplot(df, ggplot2::aes(prop, weight)) +
    ggplot2::geom_point(ggplot2::aes(colour = col)) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::xlab("Proportionality Between Features (rho)") +
    ggplot2::ylab("Discrimination Between Groups") +
    ggplot2::ggtitle("Group Discrimination by Feature Pair") +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(0, max(df$weight)) +
    ggplot2::geom_hline(yintercept = -log(.05 / nrow(df)), color = "lightgrey") +
    ggplot2::geom_hline(yintercept = -log(.05^2 / nrow(df)), color = "black")

  plot(g)

  if(!missing(k)){

    return(clust)
  }
}

#' Plot Check
#'
#' Performs data checks before plotting, triggering messages
#'  or errors when appropriate. For back-end use only.
#'
#' @inheritParams bucket
plotCheck <- function(rho, prompt){

  if(!class(rho) == "propr"){

    stop("Uh oh! You can only display a 'propr' object created by 'perb'.")
  }

  if(rho@matrix[1, 1] != 1){

    stop("Uh oh! You can only display a 'propr' object created by 'perb'.")
  }

  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Uh oh! This plot method depends on ggplot2! ",
         "Try running: install.packages('ggplot2')")
  }

  if(length(rho@pairs) != 0){

    message("Note that this display method displays all pairs, and not only indexed pairs.")
  }

  if(nrow(rho@matrix) > 1000 & prompt){

    message("Uh oh! A large number of features were detected (>1000).\n",
            "Are you sure you want to plot them all?\n",
            "0: Nevermind\n1: Proceed\n2: Hmm...")
    response <- readline(prompt = "Which do you choose? ")
    if(!response == 1) stop("Plot method aborted.")
  }
}

#' Build \code{propr} Table
#'
#' This function builds a table of VLR, VLS, and rho
#'  for each feature pair in a \code{propr} object. If the
#'  argument \code{k} is provided, the table will also
#'  include co-cluster membership (as matching the
#'  \code{\link{bucket}} or \code{\link{prism}} plots).
#'
#' @inheritParams bucket
#' @return Returns a \code{data.frame} of pairwise relationships.
#'  If the argument \code{k} is provided, returns a list of
#'  the \code{data.frame} of pairwise relationships and the
#'  cluster membership.
#'
#' @importFrom stats var
#' @export
slate <- function(rho, k, prompt = TRUE){

  if(suppressWarnings(!requireNamespace("fastcluster", quietly = TRUE))){
    stop("Uh oh! This plot method depends on fastcluster! ",
         "Try running: install.packages('fastcluster')")
  }

  plotCheck(rho, prompt)

  # Enforce some kind of colnames
  if(!is.null(colnames(rho@logratio))){ feat.each <- colnames(rho@logratio)
  }else{ feat.each <- as.character(1:ncol(rho@logratio))}

  # Calculate log-ratio transformed variances
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
#'
#' @export
prism <- function(rho, k, prompt = TRUE){

  df <- slate(rho, k, prompt)

  if(!missing(k)){

    clust <- df[[2]]
    df <- df[[1]]

  }else{

    df$CoCluster <- as.character(0)
  }

  g <-
    ggplot2::ggplot(df, ggplot2::aes_string(x = "VLS", y = "VLR")) +
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

  plot(g)

  if(!missing(k)){

    return(clust)
  }
}

#' Make Bokeh Plot
#'
#' Plots the individual feature variances for each feature
#'  pair in the \code{propr} object. Highly proportional
#'  pairs will aggregate near the \code{y = x} diagonal.
#'  Clusters that appear toward the top-right of the
#'  figure contain features with fixed abundance across
#'  all samples. Clusters that appear toward the
#'  bottom-left of the figure contain features with highly
#'  variable abundance across all samples.
#'
#' Providing the argument \code{k} will color feature pairs
#'  by co-cluster membership. In other words, a feature pair
#'  will receive a color if and only if both features belong
#'  to same the cluster (calculated using \code{hclust}).
#'
#' @inheritParams bucket
#' @return Returns cluster membership if \code{k} is provided.
#'
#' @export
bokeh <- function(rho, k, prompt = TRUE){

  df <- slate(rho, k, prompt)

  if(!missing(k)){

    clust <- df[[2]]
    df <- df[[1]]

  }else{

    df$CoCluster <- as.character(0)
  }

  df$VL1 <- -log(df$VL1)
  df$VL2 <- -log(df$VL2)

  g <-
    ggplot2::ggplot(df, ggplot2::aes_string(x = "VL1", y = "VL2")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "CoCluster", alpha = "rho")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::scale_alpha_continuous(limits = c(-1, 1), name = "Proportionality") +
    ggplot2::xlab("Log-fold *lr-transformed Variance[1]") +
    ggplot2::ylab("Log-fold *lr-transformed Variance[2]") +
    ggplot2::ggtitle("Distribution of *lr-transformed Variance") +
    ggplot2::xlim(min(df$VL1), max(df$VL1)) +
    ggplot2::ylim(min(df$VL2), max(df$VL2)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "lightgrey")

  plot(g)

  if(!missing(k)){

    return(clust)
  }
}

#' Make MDS Plot
#'
#' Plots the first two principal components as calculated
#'  using the log-ratio transformed feature vectors. This
#'  provides a statistically valid alternative to
#'  conventional principal components analysis (PCA) and
#'  multi-dimensional scaling (MDS) plotting.
#'
#' For more information, see DOI:10.1139/cjm-2015-0821.
#'
#' @inheritParams bucket
#'
#' @importFrom stats prcomp
#' @export
mds <- function(rho, group, prompt = TRUE){

  plotCheck(rho, prompt)

  if(missing(group)){

    group <- rep("None", nrow(rho@logratio))
  }

  df <- data.frame("id" = rownames(rho@logratio), "group" = as.character(group),
                   prcomp(rho@logratio)$x[, c(1, 2)])
  g <-
    ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes_string(x = "PC1", y = "PC2", colour = "group"), data = df) +
    ggplot2::geom_text(
      ggplot2::aes_string(x = "PC1", y = "PC2", label = "id", colour = "group"),
      data = df, size = 3, vjust = -1) +
    ggplot2::theme_bw() +
    ggplot2::xlab("First *lr-transformed Component") +
    ggplot2::ylab("Second *lr-transformed Component") +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Group") +
    ggplot2::ggtitle("*lr-transformed MDS Plot")

  plot(g)
}

#' Make Snapshot Plot
#'
#' Plots the log-ratio transformed feature abundance as
#'  a heatmap, with the axes ordered by dendrogram. Heatmap
#'  intensity is scaled across the feature vector to make it
#'  easier to assess visually the difference in the
#'  log-ratio transformed abundance between samples.
#'
#' @inheritParams bucket
#'
#' @importFrom stats heatmap
#' @export
snapshot <- function(rho, prompt = TRUE){

  plotCheck(rho, prompt)
  heatmap(rho@logratio, scale = "col", labCol = "")
}
