#' @rdname propd
#' @section Methods (by generic):
#' \code{plot:} Method to plot \code{propd} object.
#' @export
setMethod("plot", signature(x = "propd", y = "missing"),
          function(x, y, cutoff = 1000, col1, col2, propr, prompt = TRUE, d3 = FALSE, plotSkip = FALSE){

            names <- colnames(x@counts)
            group <- unique(x@group)

            # Apply cutoff of [0, 1] or a large integer for top theta
            if(cutoff > nrow(x@theta)) cutoff <- nrow(x@theta)
            if(cutoff > 1) cutoff <- sort(x@theta$theta)[cutoff]
            type <- x@active
            x <- x@theta[x@theta$theta <= cutoff, ]
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

            if(type == "theta_d" | type == "theta_f"){

              g <- migraph.color(g, partners[x$lrm1 > x$lrm2], pairs[x$lrm1 > x$lrm2], "coral1") # red
              g <- migraph.color(g, partners[x$lrm1 < x$lrm2], pairs[x$lrm1 < x$lrm2], "lightseagreen") # blue
              message("Red: Pair has higher LRM in group ", group[1], " than in group ", group[2])
              message("Blue: Pair has lower LRM in group ", group[1], " than in group ", group[2])

            }else if(type == "theta_e"){

              g <- migraph.color(g, partners[x$lrv1 < x$lrv2], pairs[x$lrv1 < x$lrv2], "coral1") # red
              g <- migraph.color(g, partners[x$lrv1 > x$lrv2], pairs[x$lrv1 > x$lrv2], "lightseagreen") # blue
              message("Red: Pair has nearly absent LRV in group ", group[1], " than in group ", group[2])
              message("Blue: Pair has nearly total LRV in group ", group[1], " than in group ", group[2])

            }else{

              stop("Provided theta type not supported.")
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

#' Calculate PALs for Pairs
#'
#' This function finds the Popular Adjacent Ligand or Self (PALs)
#'  for each feature in the \code{@@theta} slot of a \code{propd}
#'  object. Specifically, we define PALs as the adjacent node with
#'  the highest amount of connectivity. If node itself has more
#'  connectivity than any of its neighbors, it is its own PAL.
#'
#' @inheritParams propd
#' @return A named vector of PALs, ordered by decreasing
#'  connectivity of the input nodes. The names of the result
#'  refer to the input nodes themselves.
#'
#' @export
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

#' @rdname propd
#' @export
shale <- function(object, cutoff = 1000, k, prompt = TRUE, clean = FALSE){

  if(class(object) != "propd") stop("This function requires a 'propd' object.")
  if(!"lrm1" %in% colnames(object@theta)) stop("This function requires an 'lrm1' column!")
  if(!"lrm2" %in% colnames(object@theta)) stop("This function requires an 'lrm2' column!")

  # Apply cutoff of [0, 1] or a large integer for top theta
  if(cutoff > nrow(object@theta)) cutoff <- nrow(object@theta)
  if(cutoff > 1) cutoff <- sort(object@theta$theta)[cutoff]
  object@theta <- object@theta[object@theta$theta <= cutoff, ]
  df <- object@theta
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

#' @rdname propd
#' @export
geyser <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

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
    ggplot2::geom_abline(slope = 1/8, intercept = 0, color = "lightgrey")

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))

  }else{

    plot(g)
  }

  return(g)
}

#' @rdname propd
#' @export
bowtie <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

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
    ggplot2::geom_abline(slope = 1/8, intercept = 0, color = "lightgrey")

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))

  }else{

    plot(g)
  }

  return(g)
}

#' @rdname propd
#' @export
gemini <- function(object, cutoff = 1000, k = 5, prompt = TRUE, plotly = FALSE){

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
    ggplot2::geom_vline(xintercept = 0, color = "lightgrey")

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))

  }else{

    plot(g)
  }

  return(g)
}

#' @rdname propd
#' @export
slice <- function(object, cutoff = 1000, reference, prompt = TRUE, plotly = FALSE){

  if(missing(reference)) stop("Please provide a reference feature by name.")
  if(!is.character(reference)) stop("Please provide a reference feature by name.")
  df <- shale(object, cutoff, prompt = prompt)

  # Retrieve log-ratio counts for each pair containing the reference
  g1 <- object@group == unique(object@group)[1]
  all <-
    lapply(1:nrow(df),
           function(i){

             if(any(reference == df[i, c("Pair", "PartnerName")])){

               if(df$PairName[i] == reference){

                 bot <- df$PairName[i]
                 top <- df$PartnerName[i]

               }else if(df$PartnerName[i] == reference){

                 bot <- df$PartnerName[i]
                 top <- df$PairName[i]
               }

               lr1 <- log(object@counts[g1, top] / object@counts[g1, bot])
               lr2 <- log(object@counts[!g1, top] / object@counts[!g1, bot])

               data.frame(
                 top, bot,
                 "Group" = c(object@group[g1], object@group[!g1]),
                 "lr" = c(lr1, lr2),
                 "lrm" = c(mean(lr1), mean(lr2)),
                 "lrv" = c(var(lr1), var(lr2))
               )
             }
           })

  # Clean data for ggplot2
  final <- do.call("rbind", all)
  final$label <- colnames(object@counts)[final$top]

  g <-
    ggplot2::ggplot(final, ggplot2::aes_string(x = "label", y = "lr", LRM = "lrm", LRV = "lrv")) +
    ggplot2::geom_point(ggplot2::aes_string(col = "Group")) +
    ggplot2::theme_bw() + ggplot2::scale_colour_brewer(palette = "Set2", name = "Group") +
    ggplot2::xlab("Feature[i]") + ggplot2::ylab("Log-Ratio Abundance ( Feature[i] / Reference )") +
    ggplot2::ggtitle(paste("Log-Ratio Abundances for Reference Slice:", reference))

  if(plotly){

    packageCheck("plotly")
    return(plotly::ggplotly(g))

  }else{

    plot(g)
  }

  return(g)
}

#' @rdname propd
#' @export
decomposed <- function(object, cutoff = 1000, prompt = TRUE){

  packageCheck("compositions")
  df <- shale(object, cutoff = cutoff, prompt = prompt, clean = TRUE)

  # Generalized decomposition of LRV for weighted theta types
  decomp <- matrix(0, nrow = nrow(df), ncol = 3)
  decomp[, 1] <- df$p1 * df$LRV1 / (df$p * df$LRV)
  decomp[, 2] <- df$p2 * df$LRV2 / (df$p * df$LRV)
  decomp[, 3] <- df$p1 * df$p2 * (df$LRM2 - df$LRM1)^2 / (df$p^2 * df$LRV)

  x <- compositions::acomp(decomp)
  plot(x, pch = 20, col = grDevices::rgb(0.1, 0.1, 0.1, 0.1),
       labels = c("group 1  ", "  group 2", "between-group"), axes = TRUE)

  return(TRUE)
}
