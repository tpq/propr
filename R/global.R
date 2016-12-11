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

  }else{

    if(class(rho) != "propr"){
      stop("Uh oh! This function requires a 'propr' object.")
    }
  }

  if(is.null(colnames(rho@logratio))){
    colnames(rho@logratio) <- as.character(1:ncol(rho@logratio))
  }

  if(is.null(rownames(rho@logratio))){
    rownames(rho@logratio) <- as.character(1:nrow(rho@logratio))
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

#' Differential Proportionality Check
#'
#' Performs data checks when comparing two \code{propr} objects,
#'  triggering messages or errors when appropriate. For back-end
#'  use only.
#'
#' @inheritParams prop2prob
#' @param forceBoth Toggles whether to perform checks for the
#'  second \code{propr} object.
differentialCheck <- function(x, y, forceBoth){

  packageCheck("data.table")
  x <- suppressMessages(plotCheck(x, prompt = FALSE, plotly = FALSE, indexNaive = TRUE))

  if(!missing(y) | forceBoth){

    y <- suppressMessages(plotCheck(y, prompt = FALSE, plotly = FALSE, indexNaive = TRUE))
    if(!identical(colnames(x@logratio), colnames(y@logratio))){
      stop("Uh oh! Make sure both 'propr' objects have the same features.")
    }
  }
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
#' @source <DOI:10.1111/mec.13184>
#' @usage data(caneToad.counts)
"caneToad.counts"

#' Example Cane Toad Group Data
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
#' Includes cane toad transcripts with at least 10 counts
#'  across at least 10 samples. Used for vignette.
#'
#' @source <DOI:10.1111/mec.13184>
#' @usage data(top)
"top"
