#' Plot Check
#'
#' Performs data checks before plotting, triggering messages
#'  or errors when appropriate. For back-end use only.
#'
#' @param x An object of class \code{propr}.
plotCheck <- function(x){

  ###
  # Add check for rho
  ###

  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Uh oh! This plot method depends on ggplot2! ",
         "Try running: install.packages('ggplot2')")
  }

  if(!class(x) == "propr"){

    stop("Uh oh! You can only simplify an indexed 'propr' object.")
  }

  if(length(x@pairs) == 0){

    message("Note that this plot method plots all pairs, and not only indexed pairs.")
  }

  if(nrow(x@matrix) > 1000){

    message("Uh oh! A large number of features were detected (>1000). Are you sure you want to plot them all?\n",
            "0: Nevermind\n1: Proceed\n2: Hmm...")
    response <- readline(prompt="Which do you choose? ")
    if(!response == 1) stop("Plot method aborted.")
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
#' @param A \code{propr} object created from \code{\link{perb}}.
#' @param group. A character vector. Group or sub-group memberships,
#'  ordered according to the column names in \code{@@counts} and
#'  \code{@@logratio}. Parameter required for \code{\link{bucket}}
#'  and optional for \code{\link{mds}}.
#' @param k A numeric scalar. The number of groups into which to
#'  cluster the subjects. Clusters calculated  using \code{hclust}.
#'  Optional parameter for \code{\link{bucket}} and
#'  \code{\link{prism}}.
#'
#' @return Returns cluster membership if \code{k} is provided.
#'
#' @importFrom stats as.formula lm aov cutree hclust dist
#' @export
bucket <- function(rho, group, k){ # pronounced bouquet

  plotCheck(rho)

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

    clust <- cutree(hclust(dist(rho@matrix)), k = k)
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
    ggplot2::ggplot(df, aes(prop, weight)) +
    ggplot2::geom_point(aes(colour = col)) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::xlab("Proportionality between Features (rho)") +
    ggplot2::ylab("Discrimination between Groups") +
    ggplot2::ggtitle("Distribution of lr*-transformed Variance") +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(0, max(df$weight)) +
    ggplot2::geom_hline(yintercept = -log(.05 / nrow(df)), color = "lightgrey") +
    ggplot2::geom_hline(yintercept = -log(.05^2 / nrow(df)), color = "black")

  plot(g)

  if(!missing(k)){

    return(clust)
  }
}

#' Make Prism Plot
#'
#' Plots the variance of ratio of the log-ratio transformed feature
#'  pair (vlr) versus the sum of the individual variances of each
#'  log-ratio transformed feature (vls). The ratio of the vlr to
#'  the vls equals \code{1 - rho}. As such, we use here seven
#'  rainbow colored lines to indicate where \code{rho} equals
#'  \code{[.01, .05, .50, 0, 1.50, 1.95, 1.99]}, going from
#'  red to violet.
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
prism <- function(rho, k){

  plotCheck(rho)

  # Calculate log-ratio transformed variances
  var.ratio <- propr:::vlrRcpp(rho@counts[])
  var.each <- apply(rho@logratio, 2, var)

  # Cluster if k is provided
  if(!missing(k)){

    clust <- cutree(hclust(dist(rho@matrix)), k = k)
  }

  # Build graph components
  llt <- ncol(var.ratio) * (ncol(var.ratio)-1) * 1/2 # size of lower-left triangle
  vlr <- vector("numeric", llt) # var log ratio
  vls <- vector("numeric", llt) # var log sum
  col <- vector("numeric", llt)
  count <- 1
  for(j in 2:nrow(var.ratio)){
    for(i in 1:(j-1)){

      vlr[count] <- var.ratio[j, i]
      vls[count] <- var.each[j] + var.each[i]

      # Since each col initializes as zero
      if(!missing(k)) if(clust[i] == clust[j]) col[count] <- clust[i]

      count <- count + 1
    }
  }

  df <- data.frame(vlr, vls, "col" = as.character(col))
  g <-
    ggplot2::ggplot(df, aes(vls, vlr)) +
    ggplot2::geom_point(aes(colour = col)) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Co-Cluster") +
    ggplot2::xlab("Variance of the Log Sum (vls)") +
    ggplot2::ylab("Variance of the Log Ratio (vlr)") +
    ggplot2::ggtitle("Distribution of lr*-transformed Variance") +
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
#' @export
mds <- function(rho, group){

  plotCheck(rho)

  if(missing(group)){

    group <- rep("None", nrow(rho@logratio))
  }

  df <- data.frame("id" = rownames(rho@logratio), "group" = as.character(group),
                   prcomp(rho@logratio)$x[, c(1, 2)])
  g <-
    ggplot2::ggplot() +
    ggplot2::geom_point(aes(x = PC1, y = PC2, colour = group), data = df) +
    ggplot2::geom_text(aes(x = PC1, y = PC2, label = id, colour = group), data = df, size = 3, vjust = -1) +
    ggplot2::theme_bw() +
    ggplot2::xlab("First *lr-transformed component") +
    ggplot2::ylab("Second *lr-transformed component") +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Group") +
    ggplot2::ggtitle("*lr-transformed MDS Plot")

  plot(g)

  return(TRUE)
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
#' @export
snapshot <- function(rho){

  plotCheck(rho)
  heatmap(rho@logratio, scale = "col", labCol = "")
}
