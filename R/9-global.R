#' Update FDR by Permutation
#'
#' This function updates FDR for a set of cutoffs.
#'
#' This function wraps \code{updateCutoffs.propr} and
#'  \code{updateCutoffs.propd}.
#'
#' @inheritParams all
#' @return A \code{propr} or \code{propd} object.
#'
#' @export
updateCutoffs <- function(object, cutoff = seq(.05, .95, .3), ncores = 1){

  if(class(object) == "propr"){

    if(ncores == 1){
      message("Alert: Try parallelizing updateCutoffs with ncores > 1.")
    }

    updateCutoffs.propr(object, cutoff, ncores)

  }else if(class(object) == "propd"){

    if(ncores > 1){
      message("Alert: Parallel updateCutoffs not yet supported.")
    }

    updateCutoffs.propd(object, cutoff)

  }else{

    stop("Provided 'object' not recognized.")
  }
}

#' Make Long Data from Wide Data
#'
#' @param wide A data set in wide format.
#'
#' @return A data set in long format.
#'
#' @export
wide2long <- function(wide){

  # Force column names
  if(is.null(colnames(wide))){
    colnames(wide) <- paste0("Col", 1:ncol(wide))
  }

  # Force row names
  if(is.null(rownames(wide))){
    rownames(wide) <- as.character(rownames(wide))
  }

  df <- data.frame("value" = as.vector(as.matrix(wide)))
  df$variable <- unlist(lapply(colnames(wide), function(x) rep(x, nrow(wide))))
  df$id <- rownames(wide)
  return(df)
}

#' Recast Matrix as Feature (Log-)Ratios
#'
#' The \code{ratios} function recasts a matrix with N feature columns
#'  as a new matrix with N * (N - 1) / 2 feature (log-)ratio columns.
#'
#' When the \code{alpha} argument is provided, this function returns
#'  the (log-)ratios as \code{(partner^alpha - pair^alpha) / alpha}.
#'
#' @param matrix A matrix. The data to recast.
#' @param alpha A double. See vignette for details. Leave missing
#'  to skip Box-Cox transformation.
#' @return A matrix of (log-)ratios.
#' @export
ratios <- function(matrix, alpha = NA){

  lab <- labRcpp(ncol(matrix))

  # Replace count zeros with 1 if appropriate
  if(any(as.matrix(matrix) == 0) & is.na(alpha)){
    message("Alert: Replacing 0s with next smallest value.")
    zeros <- matrix == 0
    matrix[zeros] <- min(matrix[!zeros])
  }

  # Get (log-)ratios [based on alpha]
  if(is.na(alpha)){
    ratios <- log(matrix[, lab$Partner] / matrix[, lab$Pair])
  }else{
    message("Alert: Using alpha transformation to approximate log-ratios.")
    ratios <- (matrix[, lab$Partner]^alpha - matrix[, lab$Pair]^alpha) / alpha
  }

  # Name columns
  if(!is.null(colnames(matrix))){
    colnames(ratios) <-
      paste0(colnames(matrix)[lab$Partner],
             "/", colnames(matrix)[lab$Pair])
  }

  return(ratios)
}

#' \code{igraph} Helper Functions
#'
#' @param g An \code{igraph} object.
#' @param names1,names2 A character vector. The \code{names1}
#'  argument defines a first set of vertices. The \code{names2}
#'  argument defines a second set of vertices to which the
#'  first set connects (i.e., element-wise), thereby defining
#'  a set of edges.
#' @param col A character string. The color applied to all
#'  vertices (or edges) specified by the \code{names1} (or
#'  \code{names2}) argument.
#' @param force A boolean. If true, the function adds any
#'  missing vertices before adding edges. If false, the
#'  function only adds edges that have both vertices
#'  already present.
#'
#' @return An \code{igraph} object.
#' @name migraph

#' @rdname migraph
migraph.add <- function(g, names1, names2, force = TRUE){

  packageCheck("igraph")

  if(missing(names2)){ names2 <- NULL
  }else{ names2 <- as.character(names2) }
  names1 <- as.character(names1)

  if(force | is.null(names2)){

    # Add any missing vertices before adding edges
    all <- union(names1, names2)
    new <- all[!all %in% igraph::V(g)$name]
    if(length(new) > 0){
      g <- igraph::add.vertices(g, length(new), "name" = new)
    }

  }else{

    # Only add edges in which both vertices appear on graph
    keep <- names1 %in% igraph::V(g)$name & names2 %in% igraph::V(g)$name
    names1 <- names1[keep]
    names2 <- names2[keep]
    if(length(names1) == 0){
      stop("No new edges to add.")
    }
  }

  if(!is.null(names2)){

    if(length(names1) != length(names2)){
      stop("Argument 'names1' and 'names2' should have same length.")
    }

    edges <- unlist(
      lapply(1:length(names1),
             function(i) c(names1[i], names2[i])))

    g <- igraph::add.edges(g, edges, color = "black")
  }

  g <- igraph::as.undirected(g)
  g <- igraph::simplify(g)

  return(g)
}

#' @rdname migraph
migraph.color <- function(g, names1, names2, col){

  packageCheck("igraph")

  if(is.null(igraph::V(g)$color)) igraph::V(g)$color <- "white"
  if(is.null(igraph::E(g)$color)) igraph::E(g)$color <- "black"

  if(missing(names2)){ names2 <- NULL
  }else{ names2 <- as.character(names2) }
  names1 <- as.character(names1)

  if(is.null(names2)){

    igraph::V(g)$color <- ifelse(igraph::V(g)$name %in% names1,
                                 col, igraph::V(g)$color)

  }else{

    if(length(names1) != length(names2)){
      stop("Argument 'names1' and 'names2' should have same length.")
    }

    names <- unlist(
      lapply(1:length(names1),
             function(i) paste(names1[i], names2[i], sep = "-")))
    also <- unlist(
      lapply(1:length(names1),
             function(i) paste(names2[i], names1[i], sep = "-")))

    edges <- apply(igraph::get.edgelist(g), 1, paste, collapse = "-")
    igraph::E(g)$color <- ifelse(edges %in% names | edges %in% also,
                                 col, igraph::E(g)$color)
  }

  return(g)
}

#' @rdname migraph
migraph.clean <- function(g){

  packageCheck("igraph")

  igraph::V(g)$size <- 2
  igraph::V(g)$label <- NA

  return(g)
}

#' Make Progress Bar
#'
#' @param i The current iteration.
#' @param k Total iterations.
#' @param numTicks The result of \code{progress}.
#' @return The next \code{numTicks} argument.
progress <- function(i, k, numTicks){

  if(i == 1) numTicks <- 0

  if(numTicks == 0) cat("|-")

  while(i > numTicks*(k/40)){

    cat("-")
    if(numTicks == 10) cat("(25%)")
    if(numTicks == 20) cat("(50%)")
    if(numTicks == 30) cat("(75%)")
    numTicks <- numTicks + 1
  }

  if(i == k) cat("-|\n")

  return(numTicks)
}

#' Plot Multiple Graphs
#'
#' Easily plot multiple graphs within the same window. Code adapted from
#'  http://www.cookbook-r.com/. For back-end use only.
#'
#' @param ... Multiple plots.
#' @param cols A numeric scalar. The number of plot columns.
multiplot <- function(..., cols = 1){

  # Layout a list of plots
  plots <- list(...)
  numPlots <- length(plots)
  layout <- matrix(seq(from = 1, to = cols * ceiling(numPlots/cols)),
                   ncol = cols, nrow = ceiling(numPlots/cols))

  # Set up the page
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

  # Place each plot, in the correct location
  for(i in 1:numPlots){

    # Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    print(plots[[i]], vp = grid::viewport(
      layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
  }
}

#' Dendrogram Plot Wrapper
#'
#' Builds \code{ggplot2} dendrograms. For back-end use only.
#'
#' @param dendrogram A result from \code{as.dendrogram}.
#' @return Returns a \code{ggplot} object.
ggdend <- function(dendrogram){

  dg <- ggdendro::dendro_data(dendrogram)
  df <- dg$segments

  g <-
    ggplot2::ggplot() + ggplot2::geom_segment(
      ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), data = df) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(axis.text = ggplot2::element_blank())

  return(g)
}

#' Plot Check
#'
#' Performs data checks before plotting, triggering messages
#'  or errors when appropriate. For back-end use only.
#'
#' @inheritParams slate
#' @param indexNaive Toggles whether to perform checks for an
#'  "index-naive" plot function.
#' @return Returns a \code{propr} object with guaranteed
#'  column names and row names.
plotCheck <- function(rho, prompt, plotly, indexNaive){

  if(plotly) packageCheck("plotly")

  if(indexNaive){

    if(class(rho) != "propr" | rho@matrix[1, 1] != 1){
      stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
    }

    if(length(rho@pairs) != 0) message("Alert: This function ignores index.")
    if(prompt) promptCheck(nrow(rho@matrix))
  }

  if(class(rho) != "propr"){
    stop("Uh oh! This function requires a 'propr' object.")
  }

  return(rho)
}

#' Dendrogram Plot Check
#'
#' Performs data checks before plotting, triggering messages
#'  or errors when appropriate. For back-end use only.
dendroCheck <- function(){

  packageCheck("reshape2")
  packageCheck("ggdendro")
  packageCheck("grid")
}

#' Feature Check
#'
#' Prompts user when performing an operation on an unusually
#'  large set of features. For back-end use only.
#'
#' @param N An integer. The number of features.
promptCheck <- function(N){

  if(N > 1000){

    message("Uh oh! A large number of features were detected (>1000).\n",
            "Are you sure you want to plot them all?\n",
            "0: Nevermind\n1: Proceed\n2: Hmm...")

    response <- readline(prompt = "Which do you choose? ")
    if(!response == 1) stop("Plot method aborted.")
  }
}

#' Package Check
#'
#' Checks whether the user has the required package installed.
#'  For back-end use only.
#'
#' @param package A character string. An R package.
packageCheck <- function(package){

  if(!requireNamespace(package, quietly = TRUE)){
    stop("Uh oh! This propr method depends on ", package, "! ",
         "Try running: install.packages('", package, "')")
  }
}

#' Example Cane Toad Count Data
#'
#' Raw RNA-seq counts for twenty toads sampled from one
#'  of two regions in Australia. Data set reduced to
#'  exclude any transcripts without at least 10 counts
#'  in at least 10 samples.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(caneToad.counts)
"caneToad.counts"

#' Example Cane Toad Group Data
#'
#' Group labels for twenty toads sampled from one
#'  of two regions in Australia. Data set reduced to
#'  exclude any transcripts without at least 10 counts
#'  in at least 10 samples.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(caneToad.groups)
"caneToad.groups"

#' Mock Mail Count Data
#'
#' Includes mock count data for 5 days and 4 zip codes.
#'
#' @usage data(mail)
"mail"

#' Example propr Object
#'
#' Includes the non-transformed count abundances from all
#'  cane toad transcripts with at least 10 counts in at least
#'  10 samples, subsetted to include only those indexed
#'  by \code{rho > .995}. Used for vignette.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(top.counts)
"top.counts"

#' Example propr Object
#'
#' Includes the log-ratio transformed abundances from all
#'  cane toad transcripts with at least 10 counts in at least
#'  10 samples, subsetted to include only those indexed
#'  by \code{rho > .995}. Used for vignette.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(top.lr)
"top.lr"

#' Example Absolute mRNA
#'
#' Data generated with supplemental script provided by
#'  <DOI:10.1371/journal.pcbi.1004075>. Data originally
#'  sourced from <DOI:10.1016/j.cell.2012.09.019>.
#'  A time series of yeast mRNA abundance after removal
#'  of a key nutrient. Absolute abundance estimated
#'  by multiplying microarray signal (relative to first
#'  time point) by the initial nCounter-calibrated and
#'  copy-per-cell-adjusted RNA-seq abundance (averaged
#'  across two replicates). Divide absolute abundances
#'  by total sample abundance to make data relative.
#'
#' @usage data(marg.abs)
"marg.abs"

#' Example propd Object
#'
#' Includes results from \code{\link{propd}} as applied to
#'  cane toad transcripts with at least 40 counts in at least
#'  20 samples (after removing any transcripts with 0 counts).
#'  The resultant object is filtered to include only the top
#'  1000 theta_d values in the \code{@@results} slot.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(pd.d)
"pd.d"

#' Example propd Object
#'
#' Includes results from \code{\link{propd}} as applied to
#'  cane toad transcripts with at least 40 counts in at least
#'  20 samples (after removing any transcripts with 0 counts).
#'  The resultant object is filtered to include only the top
#'  1000 theta_e values in the \code{@@results} slot.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(pd.e)
"pd.e"
