#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about calculating proportionality, see
#'  Details.
#'
#' To learn more about visualizing proportionality, see
#'  \code{\link{visualize}}.
#'
#' To learn more about \code{ALDEx2} package integration, see
#'  \code{\link{aldex2propr}}.
#'
#' To learn more about differential proportionality, see
#'  \code{\link{propd}}.
#'
#' To learn more about compositional data analysis, and its relevance
#'  to biological count data, see the bundled vignette.
#'
#' @details
#' Let D represent a number of features measured across N samples.
#' 	This function calculates proportionality from
#' 	a data set with N rows and D columns.
#'  One can think of phi as
#' 	analogous to a distance matrix, except that it has no symmetry unless forced.
#' 	One can think of rho as
#' 	analogous to a correlation matrix.
#' 	One can think of phs as
#' 	either a naturally symmetric variant of phi or a monotonic variant of rho.
#' 	Also, one can use \code{corr}
#' 	to calculate correlation from log-ratio transformed data.
#'
#' This function depends on a reference and uses the centered log-ratio
#'  transformation by default. The user may also specify any number of
#'  features (by index or name) to use as a reference instead.
#'  Alternatively, \code{ivar = "iqlr"} will transform data using the
#'  geometric mean of features with variances that fall in the
#'  inter-quartile range of all per-feature variances (based on
#'  the \code{ALDEx2} package).
#'
#' The \code{propr} method calculates proportionality. This fails in
#'  the setting of zero counts. The \code{propr} method
#'  will use a Box-Cox transformation to approximate VLR based on
#'  the parameter \eqn{\alpha}, if provided. We refer the user to
#'  the vignette for more details.
#'
#' @slot counts A data.frame. Stores the original "count matrix" input.
#' @slot alpha A double. Stores the alpha value used for transformation.
#' @slot metric A character string. The metric used to calculate proportionality.
#' @slot ivar A vector. The reference used to calculate proportionality.
#' @slot logratio A data.frame. Stores the transformed "count matrix".
#' @slot matrix A matrix. Stores the proportionality matrix.
#' @slot pairs A vector. Indexes the proportional pairs of interest.
#' @slot results A data.frame. Stores the pairwise \code{propr} measurements.
#' @slot permutes A list. Stores the shuffled transformed "count matrix"
#'  instances, used to reproduce permutations of \code{propr}.
#' @slot fdr A data.frame. Stores the FDR cutoffs for \code{propr}.
#'
#' @inheritParams all
#' @return Returns a \code{propr} object.
#'
#' @name propr
#' @useDynLib propr, .registration = TRUE
#' @importFrom methods show new
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname propr
#' @export
setClass("propr",
         slots = c(
           counts = "data.frame",
           alpha = "numeric",
           metric = "character",
           ivar = "ANY",
           logratio = "data.frame",
           matrix = "matrix",
           pairs = "numeric",
           results = "data.frame",
           permutes = "list",
           fdr = "data.frame"
         )
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propr} object.
#' @export
setMethod("show", "propr",
          function(object){

            cat("Not weighted", "and",
                ifelse(is.na(object@alpha), "not alpha-transformed", "alpha-transformed"), "\n")

            cat("@counts summary:",
                nrow(object@counts), "subjects by", ncol(object@counts), "features\n")

            cat("@logratio summary:",
                nrow(object@logratio), "subjects by", ncol(object@logratio), "features\n")

            cat("@matrix summary:",
                nrow(object@matrix), "features by", ncol(object@matrix), "features\n")

            if(length(object@pairs) > 0 | nrow(object@matrix) == 0){

              cat("@pairs summary:", length(object@pairs), "feature pairs\n")

            }else{

              cat("@pairs summary: index with `[` method\n")
            }

            cat("@fdr summary:",
                ncol(object@permutes), "iterations\n")

            if(nrow(object@fdr) > 0){
              print(object@fdr)
            }

            cat("See ?propr for object methods\n")
          }
)

#' @rdname propr
#' @export
propr <- function(counts, metric = c("rho", "phi", "phs", "cor", "vlr"), ivar = "clr",
                  select, symmetrize = FALSE, alpha, p = 100){

  # Clean "count matrix"
  # if(any(apply(counts, 2, function(x) all(x == 0)))){
  #   stop("Remove components with all zeros before proceeding.")}
  if(any(counts < 0)) stop("Data may not contain negative measurements.")
  if(any(is.na(counts))) stop("Remove NAs from 'counts' before proceeding.")
  if(class(counts) == "data.frame") counts <- as.matrix(counts)
  if(is.null(colnames(counts))) colnames(counts) <- as.character(1:ncol(counts))
  if(is.null(rownames(counts))) rownames(counts) <- as.character(1:nrow(counts))
  if(missing(alpha)){ alpha <- NA
  }else{ if(!is.na(alpha)) if(alpha == 0) alpha <- NA }
  ct <- counts

  # Replace zeros unless alpha is provided
  if(any(as.matrix(counts) == 0) & is.na(alpha)){
    message("Alert: Replacing 0s with next smallest value.")
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])
  }

  # Establish reference
  use <- ivar2index(ct, ivar)

  # Transform counts based on reference
  if(is.na(alpha)){

    # Use g(x) = Mean[log(x)] to log-ratio transform data
    message("Alert: Saving log-ratio transformed counts to @logratio.")
    logX <- log(ct)
    logSet <- logX[, use, drop = FALSE]
    ref <- rowMeans(logSet)
    lr <- sweep(logX, 1, ref, "-")

  }else{

    # Since y = log(x) = [x^a-1]/a, ref = Mean[log(x)] = Mean[y]
    # Calculate log(x/ref) = log(x) - log(ref) = [x^a-1]/a - [ref^a-1]/a
    message("Alert: Saving alpha transformed counts to @logratio.")
    aX <- (ct^alpha - 1) / alpha
    aSet <- aX[, use, drop = FALSE]
    ref <- rowMeans(aSet)
    lr <- sweep(aX, 1, ref, "-")
  }

  # Optionally apply select for performance gain
  if(!missing(select)){

    message("Alert: Using 'select' disables permutation testing.")
    p <- 0

    # Make select boolean (it's OK if it's integer)
    if(!is.vector(select)) stop("Provide 'select' as vector.")
    if(is.character(select)) select <- match(select, colnames(counts))
    if(any(is.na(select))) stop("Some 'select' not in data.")

    counts <- counts[, select]
    ct <- ct[, select]
    lr <- lr[, select]
  }

  # Calculate proportionality
  lrv <- lr2vlr(lr)
  metric <- metric[1]
  if(metric == "rho"){
    mat <- lr2rho(lr)
  }else if(metric == "phi"){
    mat <- lr2phi(lr)
    if(symmetrize) symRcpp(mat) # optionally force symmetry
  }else if(metric == "phs"){
    mat <- lr2phs(lr)
  }else if(metric == "cor"){
    mat <- stats::cor(lr)
  }else if(metric == "vlr"){
    mat <- lrv
  }else{
    stop("Provided 'metric' not recognized.")
  }

  # Build propr object
  result <- new("propr")
  result@counts <- as.data.frame(ct)
  if(!missing(alpha)){ result@alpha <- as.numeric(alpha)
  }else{ result@alpha <- as.numeric(NA) }
  result@metric <- metric[1]
  result@ivar <- ivar
  result@logratio <- as.data.frame(lr)
  result@pairs <- vector("numeric")

  # Clean row and column names
  result@matrix <- mat
  colnames(result@matrix) <- colnames(result@logratio)
  rownames(result@matrix) <- colnames(result@logratio)

  # Set up @permutes
  result@permutes <- list(NULL)
  if(p > 0){
    # Sample compositions so that the clr does not change
    message("Alert: Fixing permutations to active random seed.")
    permutes <- vector("list", p)
    for(ins in 1:length(permutes)) permutes[[ins]] <- t(apply(ct, 1, sample))
    result@permutes <- permutes
  }

  # Set up @results
  labels <- labRcpp(ncol(lr))
  result@results <-
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      "lrv" = lltRcpp(lrv),
      "metric" = factor(metric),
      "alpha" = factor(alpha),
      "propr" = lltRcpp(mat)
    )

  # Initialize @results -- Tally frequency of 0 counts
  if(any(as.matrix(counts) == 0)){
    message("Alert: Tabulating the presence of 0 counts.")
    result@results$Zeros <- ctzRcpp(as.matrix(counts)) # count 0s
  }

  message("Alert: Use '[' to index proportionality matrix.")
  message("Alert: Use 'updateCutoffs' to calculate FDR.")

  return(result)
}

#' @rdname propr
#' @section Functions:
#' \code{phit:}
#'  A wrapper for \code{propr(counts, metric = "phi", ...)}.
#' @export
phit <- function(counts, ...){
  propr(counts, metric = "phi", ...)
}

#' @rdname propr
#' @section Functions:
#' \code{perb:}
#'  A wrapper for \code{propr(counts, metric = "rho", ...)}.
#' @export
perb <- function(counts, ...){
  propr(counts, metric = "rho", ...)
}

#' @rdname propr
#' @section Functions:
#' \code{phis:}
#'  A wrapper for \code{propr(counts, metric = "phs", ...)}.
#' @export
phis <- function(counts, ...){
  propr(counts, metric = "phs", ...)
}

#' @rdname propr
#' @section Functions:
#' \code{corr:}
#'  A wrapper for \code{propr(counts, metric = "cor", ...)}.
#' @export
corr <- function(counts, ...){
  propr(counts, metric = "cor", ...)
}
