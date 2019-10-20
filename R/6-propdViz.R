#' Calculate PALs for Pairs
#'
#' This function finds the Popular Adjacent Ligand or Self (PALs)
#'  for each feature in the \code{@@results} slot of a \code{propd}
#'  object. Specifically, we define PALs as the adjacent node with
#'  the highest amount of connectivity. If node itself has more
#'  connectivity than any of its neighbors, it is its own PAL.
#'
#' @inheritParams all
#' @return A named vector of PALs, ordered by decreasing
#'  connectivity of the input nodes. The names of the result
#'  refer to the input nodes themselves.
pals <- function(object, k){

  # Get entire graph object
  g <- suppressMessages(plot(object, cutoff = 1, prompt = FALSE, plotSkip = TRUE))

  # Get graph parameters
  d <- sort(igraph::degree(g), decreasing = TRUE)
  el <- igraph::get.edgelist(g)

  pal <- vector("integer", length(d))
  names(pal) <- names(d)
  for(i in 1:length(d)){

    # Extract all partners for the i-th object
    node.i <- names(d)[i]
    el.i <- el[node.i == el[, 1] | node.i == el[, 2], , drop = FALSE]
    z <- union(el.i[, 1], el.i[, 2])

    # Index most popular partner
    pal[i] <- names(d)[names(d) %in% z][1]
  }

  if(!missing(k)){

    # Only keep the k largest PALs
    keep <- names(sort(table(pal), decreasing = TRUE))[1:k]
    pal[!pal %in% keep] <- "Missing"
  }

  return(pal)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{propd2propr:}
#'  Transforms a \code{propd} object into a \code{propr} object
#'  where the \code{@@matrix} slot contains \eqn{1 - \theta}.
#'  Allows the user to interrogate theta using any
#'  visualization built for \code{propr} objects.
#' @rdname propd
#' @export
propd2propr <- function(object, ivar){

  prop <- new("propr")
  prop@matrix <- 1 - half2mat(object@results$theta)
  return(prop)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{plot:}
#'  Plots the interactions between pairs as a network.
#'  When plotting disjointed proportionality, red edges
#'  indicate that LRM1 > LRM2 while blue edges indicate
#'  that LRM1 < LRM2. When plotting emergent proportionality,
#'  red edges indicate that VLR1 < VLR2 while blue edges
#'  indicate that VLR1 > VLR2. Group labels numbered based on
#'  the order of the \code{group} argument to \code{propd}.
#'  Use \code{col1} and \code{col2} arguments to color nodes.
#'  For more control over the visualization of the network,
#'  consider exporting the table from \code{shale} to a
#'  network visualization tool like Cytoscape.
#' @export
setMethod("plot", signature(x = "propd", y = "missing"),
          function(x, y, cutoff = 1000, col1, col2, propr, prompt = TRUE, d3 = FALSE, plotSkip = FALSE){

            if(length(unique(x@group)) > 2){
              stop("This method only supports the analysis of 2 groups.")
            }

            names <- colnames(x@counts)
            group <- unique(x@group)

            # Apply cutoff of [0, 1] or a large integer for top theta
            if(cutoff > nrow(x@results)) cutoff <- nrow(x@results)
            if(cutoff > 1) cutoff <- sort(x@results$theta)[cutoff]
            type <- x@active
            x <- x@results[x@results$theta <= cutoff, ]
            if(nrow(x) == 0) stop("No theta remain after cutoff.")
            if(prompt) promptCheck(nrow(x))

            partners <- names[x$Partner]
            pairs <- names[x$Pair]

            # Color edges based on direction of log-fold change in VLR
            g <- igraph::make_empty_graph(directed = FALSE)
            g <- migraph.add(g, partners, pairs)

            if(!missing(propr)){

              ff <- tempfile()
              grDevices::png(filename = ff)
              cyt <- suppressMessages(propr::cytescape(propr, prompt = FALSE))
              grDevices::dev.off()
              unlink(ff)

              g <- migraph.add(g, cyt$Partner, cyt$Pair)
              g <- migraph.color(g, cyt$Partner, cyt$Pair, "forestgreen")
              message("Green: Pair positively proportional across all samples.")

              invProp <- cyt[, 3] < 0
              if(any(invProp)){

                g <- migraph.color(g, cyt$Partner[invProp], cyt$Pair[invProp], "burlywood4")
                message("Brown: Pair inversely proportional across all samples.")
              }
            }

            if(type == "theta_e" | type == "theta_h"){

              g <- migraph.color(g, partners[x$lrv1 < x$lrv2], pairs[x$lrv1 < x$lrv2], "coral1") # red
              g <- migraph.color(g, partners[x$lrv1 > x$lrv2], pairs[x$lrv1 > x$lrv2], "lightseagreen") # blue
              message("Red: Nearly all of total LRV explained by ", group[2])
              message("Blue: Nearly all of total LRV explained by ", group[1])

            }else{

              g <- migraph.color(g, partners[x$lrm1 > x$lrm2], pairs[x$lrm1 > x$lrm2], "coral1") # red
              g <- migraph.color(g, partners[x$lrm1 < x$lrm2], pairs[x$lrm1 < x$lrm2], "lightseagreen") # blue
              message("Red: Pair has higher LRM in group ", group[1], " than in group ", group[2])
              message("Blue: Pair has higher LRM in group ", group[2], " than in group ", group[1])
            }

            # Optional coloring of nodes
            if(!missing(col1)) g <- migraph.color(g, col1, col = "darkred")
            if(!missing(col2)) g <- migraph.color(g, col2, col = "darkslateblue")
            g <- migraph.clean(g)

            # Used by pals function
            if(plotSkip) return(g)

            if(d3){

              packageCheck("rgl")
              coords <- igraph::layout_with_fr(g, dim = 3)
              suppressWarnings(igraph::rglplot(g, layout = coords))

            }else{

              plot(g)
            }

            return(g)
          }
)

#' Build \code{propd} Results Table
#'
#' Builds a table of within-group and total log-ratio
#'  variances, log-ratio means, and PALs (see: \code{\link{pals}}).
#'  If the argument \code{k} is provided, the table will
#'  label at most \code{k} top PALs. Just as each node
#'  gets assigned a PAL, \code{shale} aims to assign
#'  each edge a PAL. Edges that have a top PAL as one
#'  and only one of their nodes get assigned that PAL.
#'  Edges that have top PALs as both of their nodes
#'  get assigned "Bridged". Edges without a top PAL
#'  as one of their nodes will get assigned a PAL if
#'  either (a) both nodes have the same neighbor PAL
#'  or (b) one node has a "Missing" neighbor PAL.
#'  The \code{cutoff} argument guides the maximum value of
#'  theta above which to exclude the pair. A large integer
#'  \code{cutoff} will instead retrieve the top N pairs as
#'  ranked by theta.
#'
#' @inheritParams all
shale <- function(object, cutoff = 1000, k, prompt = TRUE, clean = FALSE){

  if(class(object) != "propd") stop("This function requires a 'propd' object.")
  if(!"lrm1" %in% colnames(object@results)) stop("This function requires an 'lrm1' column!")
  if(!"lrm2" %in% colnames(object@results)) stop("This function requires an 'lrm2' column!")

  # Apply cutoff of [0, 1] or a large integer for top theta
  if(cutoff > nrow(object@results)) cutoff <- nrow(object@results)
  if(cutoff > 1) cutoff <- sort(object@results$theta)[cutoff]
  object@results <- object@results[object@results$theta <= cutoff, ]
  df <- object@results
  if(nrow(df) == 0) stop("No theta remain after cutoff.")
  if(prompt) promptCheck(nrow(df))

  # Retrieve useful information
  colnames(df)[colnames(df) %in% c("lrv", "lrv1", "lrv2", "lrm1", "lrm2")] <-
    c("LRV", "LRV1", "LRV2", "LRM1", "LRM2")
  df$FoldLRV <- log(df$LRV1 / df$LRV2)
  df$FoldLRV[is.na(df$FoldLRV)] <- 0
  df$DiffLRM <- df$LRM1 - df$LRM2
  df$PartnerName <- colnames(object@counts)[df$Partner]
  df$PairName <- colnames(object@counts)[df$Pair]

  # Get PALs for each node
  pal <- pals(object, k)
  df$PartnerPAL <- pal[df$PartnerName]
  df$PairPAL <- pal[df$PairName]

  # Assign an edge a PAL if only one node is a top PAL
  topPALs <- unique(pal)
  df$PAL <- "Missing"
  df$PAL <- ifelse(df$PartnerName %in% topPALs, df$PartnerName, df$PAL)
  df$PAL <- ifelse(df$PairName %in% topPALs, df$PairName, df$PAL)

  # Assign "Bridged" if both nodes are a top PAL
  df$PAL <- ifelse(df$PartnerName %in% topPALs & df$PairName %in% topPALs,
                   "Bridged", df$PAL)

  # Assign PAL to edges still without a top PAL
  df$PAL <- ifelse(df$PAL == "Missing" &
                     df$PartnerPAL != "Missing" & df$PairPAL == "Missing",
                   df$PartnerPAL, df$PAL)
  df$PAL <- ifelse(df$PAL == "Missing" &
                     df$PartnerPAL == "Missing" & df$PairPAL != "Missing",
                   df$PairPAL, df$PAL)
  df$PAL <- ifelse(df$PAL == "Missing" &
                     df$PartnerPAL == df$PairPAL,
                   df$PartnerPAL, df$PAL)

  # Calculate difference in LRM relative to PAL
  df$LRM1byPAL <- ifelse(df$PairName == df$PAL,
                         df$LRM1, ifelse(df$PartnerName == df$PAL,
                                         -1 * df$LRM1, 0))
  df$LRM2byPAL <- ifelse(df$PairName == df$PAL,
                         df$LRM2, ifelse(df$PartnerName == df$PAL,
                                         -1 * df$LRM2, 0))
  df$DiffLRMbyPAL <- df$LRM1byPAL - df$LRM2byPAL

  # Clean data for plot -- used by 'geyser', 'bowtie', and 'gemini'
  if(clean) return(df[df$PairName == df$PAL | df$PartnerName == df$PAL, ])

  return(df)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{geyser:}
#'  Plots indexed pairs based on the within-group
#'  log-ratio variance (VLR) for each group. Pairs near the
#'  origin show a highly proportional relationship in
#'  both groups. Each line away from the \code{y = x} line
#'  indicates a doubling of VLR compared to the other group.
#'  All pairs colored based on PAL
#'  (see: \code{\link{pals}}).
#'  See \code{gemini}.
#' @export
geyser <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

  if(length(unique(object@group)) > 2){
    stop("This method only supports the analysis of 2 groups.")
  }

  df <- shale(object, cutoff, k, prompt, clean = TRUE)

  maxLim <- max(c(df$LRV1, df$LRV2))
  g <-
    ggplot2::ggplot(
      df, ggplot2::aes_string(Partner = "PartnerName",
                              Pair = "PairName", PartnerPAL = "PartnerPAL",
                              PairPAL = "PairPAL", theta = "theta")) +
    ggplot2::geom_point(ggplot2::aes_string(x = "LRV1", y = "LRV2", colour = "PAL")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "PAL") +
    ggplot2::xlab(paste("Log-Ratio Variance (VLR) of Pair in",
                        unique(object@group)[1])) +
    ggplot2::ylab(paste("Log-Ratio Variance (VLR) of Pair in",
                        unique(object@group)[2])) +
    ggplot2::ggtitle("Distribution of Pairs by Within-Group VLR") +
    ggplot2::xlim(0, maxLim) +
    ggplot2::ylim(0, maxLim) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 2/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/2, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 4/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/4, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 8/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/8, intercept = 0, color = "lightgrey") +
    ggplot2::theme(text = ggplot2::element_text(size=18),
                   plot.title = ggplot2::element_text(size=24))

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))
  }

  return(g)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{bowtie:}
#'  Plots indexed pairs based on the log-ratio means
#'  (LRM), relative to its PAL, for each group. Pairs near
#'  the origin show comparable LRM, relative to its PAL, in
#'  both groups. Each line away from the \code{y = x} line
#'  indicates a doubling of LRM compared to the other group.
#'  All pairs colored based on PAL
#'  (see: \code{\link{pals}}).
#'  See \code{gemini}.
#' @export
bowtie <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

  if(length(unique(object@group)) > 2){
    stop("This method only supports the analysis of 2 groups.")
  }

  df <- shale(object, cutoff, k, prompt, clean = TRUE)

  maxLim <- max(abs(c(df$LRM1byPAL, df$LRM2byPAL)))
  g <-
    ggplot2::ggplot(
      df, ggplot2::aes_string(Partner = "PartnerName",
                              Pair = "PairName", PartnerPAL = "PartnerPAL",
                              PairPAL = "PairPAL", theta = "theta")) +
    ggplot2::geom_point(ggplot2::aes_string(x = "LRM1byPAL", y = "LRM2byPAL", colour = "PAL")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "PAL") +
    ggplot2::xlab(paste("Log-Ratio Mean (LRM) Relative to PAL in",
                        unique(object@group)[1])) +
    ggplot2::ylab(paste("Log-Ratio Mean (LRM) Relative to PAL in",
                        unique(object@group)[2])) +
    ggplot2::ggtitle("Distribution of Pairs by Within-Group LRM Relative to PAL") +
    ggplot2::xlim(-1 * maxLim, maxLim) +
    ggplot2::ylim(-1 * maxLim, maxLim) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 2/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/2, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 4/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/4, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 8/1, intercept = 0, color = "lightgrey") +
    ggplot2::geom_abline(slope = 1/8, intercept = 0, color = "lightgrey") +
    ggplot2::theme(text = ggplot2::element_text(size=18),
                   plot.title = ggplot2::element_text(size=24))

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))
  }

  return(g)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{gemini:}
#'  Plots indexed pairs based on the log-fold difference
#'  in log-ratio variance (VLR) between the two groups
#'  versus the difference in log-ratio means (LRM). In this
#'  figure, the LRM for each group is signed (i.e., positive
#'  or negative) such that the PAL is the denominator
#'  of the log-ratio. This allows for a fluid comparison
#'  between pairs within the same PAL module. Pairs with a
#'  "Bridged" or "Missing" PAL get excluded from this graph.
#'  Remember that an increase in VLR suggests less
#'  proportionality. All pairs colored based on PAL
#'  (see: \code{\link{pals}}).
#' @export
gemini <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

  if(length(unique(object@group)) > 2){
    stop("This method only supports the analysis of 2 groups.")
  }

  df <- shale(object, cutoff, k, prompt, clean = TRUE)

  maxLimx <- max(abs(df$DiffLRMbyPAL))
  maxLimy <- max(abs(df$FoldLRV[df$FoldLRV != Inf & df$FoldLRV != -Inf]))
  g <-
    ggplot2::ggplot(
      df, ggplot2::aes_string(Partner = "PartnerName",
                              Pair = "PairName", PartnerPAL = "PartnerPAL",
                              PairPAL = "PairPAL", theta = "theta",
                              LRV1 = "LRV1", LRV2 = "LRV2")) +
    ggplot2::geom_point(ggplot2::aes_string(x = "DiffLRMbyPAL", y = "FoldLRV", colour = "PAL")) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "PAL") +
    ggplot2::xlab(paste0("Difference in LRM Relative to PAL (",
                         unique(object@group)[1], " - ", unique(object@group)[2], ")")) +
    ggplot2::ylab(paste0("Log-Fold Change in VLR (",
                         unique(object@group)[1], " / ", unique(object@group)[2], ")")) +
    ggplot2::ggtitle("Distribution of Pairs by VLR and LRM Difference") +
    ggplot2::xlim(-1 * maxLimx, maxLimx) +
    ggplot2::ylim(-1 * maxLimy, maxLimy) +
    ggplot2::geom_abline(slope = 0, intercept = 0, color = "lightgrey") +
    ggplot2::geom_vline(xintercept = 0, color = "lightgrey") +
    ggplot2::theme(text = ggplot2::element_text(size=18),
                   plot.title = ggplot2::element_text(size=24))

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))
  }

  return(g)
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{decomposed:}
#'  Plots the decomposition of log-ratio variance into
#'  (weighted) group variances and between-group variance.
#'  Useful for visualizing how a theta type selects pairs.
#' @export
decomposed <- function(object, cutoff = 1000){

  if(length(unique(object@group)) > 2){
    stop("This method only supports the analysis of 2 groups.")
  }

  packageCheck("compositions")
  df <- getResults(object, cutoff)

  # Generalized decomposition of LRV for weighted theta types
  decomp <- matrix(0, nrow = nrow(df), ncol = 3)
  decomp[, 1] <- df$p1 * df$lrv1 / (df$p * df$lrv)
  decomp[, 2] <- df$p2 * df$lrv2 / (df$p * df$lrv)
  decomp[, 3] <- df$p1 * df$p2 * (df$lrm2 - df$lrm1)^2 / (df$p^2 * df$lrv)

  x <- suppressWarnings(compositions::acomp(decomp))
  suppressWarnings(
    plot(x, pch = 20, col = grDevices::rgb(0.1, 0.1, 0.1, 0.1),
         labels = c("group 1  ", "  group 2", "between-group"), axes = TRUE)
  )
}

#' @rdname visualize
#' @section \code{propd} Functions:
#' \code{parallel:}
#'  Plots the sample-wise log-ratio abundance across all
#'  pairs selected by the provided cutoff. Use the
#'  \code{reference} argument to subset the plot to only
#'  include pairs that contain this reference.
#' @export
parallel <- function(object, cutoff = 1000, include = NA, or = TRUE, plotly = FALSE){

  df <- getRatios(object, cutoff, include = include, or = or, melt = TRUE)
  df$variable <- factor(df$variable, levels = unique(df$variable))
  df$group <- object@group

  g <- ggplot2::ggplot(df, ggplot2::aes_string(x = "variable", y = "value", group = "id", col = "group")) +
    ggplot2::geom_line() + ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Group") +
    ggplot2::xlab("Feature Pair") + ggplot2::ylab("Log-Ratio Abundance") +
    ggplot2::ggtitle("Sample-wise Distribution of Log-Ratios Across Pairs") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::theme(text = ggplot2::element_text(size=18),
                   plot.title = ggplot2::element_text(size=24))

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))
  }

  return(g)
}
