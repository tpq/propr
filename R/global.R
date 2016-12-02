#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about how to calculate proportionality, see
#'  see \code{\link{phit}} and \code{\link{perb}}.
#'
#' To learn more about the resultant \code{propr} class object, see
#'  see the class definition \code{\link{propr-class}}.
#'
#' To learn more about \code{propr} plots, see \code{\link{smear}},
#'  \code{\link{dendrogram}}, \code{\link{bucket}}, \code{\link{prism}},
#'  \code{\link{bokeh}}, \code{\link{mds}}, and \code{\link{snapshot}}.
#'
#' To learn more about differential proportionality, see
#'  \code{\link{prop2prob}} and \code{\link{abstract}}.
#'
#' To learn more about compositional data analysis, and its relevance
#'  to biological count data, see the bundled vignette.
#'
#' @name propr
NULL

#' An S4 class to hold results from proportionality analysis.
#'
#' @slot counts A matrix. Stores the original "count matrix" input.
#' @slot logratio A matrix. Stores the log-ratio transformed "count matrix".
#' @slot matrix A matrix. Stores the proportionality matrix calculated by
#'  \code{phiRcpp} or \code{rhoRcpp}.
#' @slot pairs A vector. Indexes the proportionality metrics of interest.
#'
#' @seealso \code{\link{propr}}, \code{\link{phit}}, \code{\link{perb}}
#'
#' @export
setClass("propr",
         slots = c(
           counts = "matrix",
           logratio = "matrix",
           matrix = "matrix",
           pairs = "numeric"
         )
)

#' @useDynLib propr
#' @importFrom Rcpp sourceCpp
NULL

#' Example Cane Toad Count Data
#' @source <DOI:10.1111/mec.13184>
#' @usage data(caneToad.counts)
"caneToad.counts"

#' Example Cane Toad Group Data
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

#' Plot Check
#'
#' Performs data checks before plotting, triggering messages
#'  or errors when appropriate. For back-end use only.
#'
#' @inheritParams bucket
#' @param indexNaive Toggles whether to perform checks for an
#'  "index-naive" plot function.
#' @return Returns a \code{propr} object with guaranteed
#'  column names and row names.
plotCheck <- function(rho, prompt, plotly, indexNaive){

  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Uh oh! This display method depends on ggplot2! ",
         "Try running: install.packages('ggplot2')")
  }

  if(plotly){

    if(!requireNamespace("plotly", quietly = TRUE)){
      stop("Uh oh! This display method depends on plotly! ",
           "Try running: install.packages('plotly')")
    }
  }

  if(!class(rho) == "propr"){

    stop("Uh oh! You can only display a 'propr' object with this function.")
  }

  if(indexNaive){

    if(!requireNamespace("fastcluster", quietly = TRUE)){
      stop("Uh oh! This display method depends on fastcluster! ",
           "Try running: install.packages('fastcluster')")
    }

    if(rho@matrix[1, 1] != 1){

      stop("Uh oh! You can only display a 'propr' object created by 'perb'.")
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

  # Force some kind of column names
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

  if(!requireNamespace("reshape2", quietly = TRUE)){
    stop("Uh oh! This display method depends on reshape2! ",
         "Try running: install.packages('reshape2')")
  }

  if(!requireNamespace("ggdendro", quietly = TRUE)){
    stop("Uh oh! This display method depends on ggdendro! ",
         "Try running: install.packages('ggdendro')")
  }

  if(!requireNamespace("fastcluster", quietly = TRUE)){
    stop("Uh oh! This display method depends on fastcluster! ",
         "Try running: install.packages('fastcluster')")
  }

  if(!requireNamespace("grid", quietly = TRUE)){
    stop("Uh oh! This display method depends on grid! ",
         "Try running: install.packages('grid')")
  }
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

  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Uh oh! This propr method depends on data.table! ",
         "Try running: install.packages('data.table')")
  }

  if(class(x) != "propr") stop("Uh oh! This function requires a 'propr' object for 'x'.")
  if(!x@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
  if(is.null(colnames(x@logratio))){
    colnames(x@logratio) <- as.character(1:ncol(x@logratio))
  }

  if(!missing(y) | forceBoth){

    if(class(y) != "propr") stop("Uh oh! This function requires a 'propr' object for 'y'.")
    if(!y@matrix[1, 1]) stop("Uh oh! This function requires a 'propr' object created by 'perb'.")
    if(is.null(colnames(y@logratio))){
      colnames(y@logratio) <- as.character(1:ncol(y@logratio))
    }

    if(!identical(colnames(x@logratio), colnames(y@logratio))){
      stop("Uh oh! Make sure both 'propr' objects have the same features.")
    }
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
    ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), data = df) +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(axis.text = ggplot2::element_blank())

  return(g)
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

  if(numPlots == 1){

    print(plots[[1]])

  }else{

    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout)))
    )

    # Place each plot, in the correct location
    for(i in 1:numPlots){

      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(
        plots[[i]],
        vp = grid::viewport(layout.pos.row = matchidx$row,
                            layout.pos.col = matchidx$col)
      )
    }
  }
}
