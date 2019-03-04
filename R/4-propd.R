#' The propd Method
#'
#' @description
#' Welcome to the \code{propd} method!
#'
#' Let \eqn{X} and \eqn{Y} be non-zero positive feature vectors
#'  measured across \eqn{N} samples belonging to one of two groups,
#'  sized \eqn{N1} and \eqn{N2}. We use VLR to denote the variance
#'  of the log of the ratio of the vectors \eqn{X} over \eqn{Y}.
#'  We define theta as the weighted sum of the within-group VLR
#'  divided by the weighted total VLR.
#'
#' The \code{propd} method calculates theta. This fails in
#'  the setting of zero counts. The \code{propd} method
#'  will use a Box-Cox transformation to approximate VLR based on
#'  the parameter \eqn{\alpha}, if provided. We refer the user to
#'  the vignette for more details.
#'
#' Note that Group 1 always refers to the first element of the
#'  \code{group} vector argument supplied to \code{propd}.
#'
#' @slot counts A data.frame. Stores the original "count matrix" input.
#' @slot alpha A double. Stores the alpha value used for transformation.
#' @slot group A character vector. Stores the original group labels.
#' @slot weighted A logical. Stores whether the theta is weighted.
#' @slot weights A matrix. If weighted, stores the limma-based weights.
#' @slot active A character. Stores the name of the active theta type.
#' @slot Fivar ANY. Stores the reference used to moderate theta.
#' @slot dfz A double. Stores the prior df used to moderate theta.
#' @slot results A data.frame. Stores the pairwise \code{propd} measurements.
#' @slot permutes A data.frame. Stores the shuffled group labels,
#'  used to reproduce permutations of \code{propd}.
#' @slot fdr A data.frame. Stores the FDR cutoffs for \code{propd}.
#'
#' @inheritParams all
#' @return Returns a \code{propr} object.
#'
#' @name propd
#' @importFrom methods show new
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname propd
#' @export
setClass("propd",
         slots = c(
           counts = "data.frame",
           alpha = "numeric",
           group = "character",
           weighted = "logical",
           weights = "matrix",
           active = "character",
           Fivar = "ANY",
           dfz = "numeric",
           results = "data.frame",
           permutes = "data.frame",
           fdr = "data.frame"
         )
)

#' @rdname propd
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propd} object.
#' @export
setMethod("show", "propd",
          function(object){

            cat(ifelse(object@weighted, "Weighted", "Not weighted"), "and",
                ifelse(is.na(object@alpha), "not alpha-transformed", "alpha-transformed"), "\n")

            cat("@counts summary:",
                nrow(object@counts), "subjects by", ncol(object@counts), "features\n")

            cat("@group summary:", length(unique(object@group)), "unique groups (",
                paste(table(object@group), collapse = " x "), ")\n")

            cat("@results summary:",
                nrow(object@results), "feature pairs (", object@active, ")\n")

            cat("@fdr summary:",
                ncol(object@permutes), "iterations\n")

            if(nrow(object@fdr) > 0){
              print(object@fdr)
            }

            cat("See ?propd for object methods\n")
          }
)

#' @rdname propd
#' @export
propd <- function(counts, group, alpha, p = 100, weighted = FALSE){

  # Clean "count matrix"
  if(any(counts < 0)) stop("Data may not contain negative measurements.")
  if(any(is.na(counts))) stop("Remove NAs from 'counts' before proceeding.")
  if(class(counts) == "data.frame") counts <- as.matrix(counts)
  if(is.null(colnames(counts))) colnames(counts) <- as.character(1:ncol(counts))
  if(is.null(rownames(counts))) rownames(counts) <- as.character(1:nrow(counts))

  # Clean group
  if(length(unique(group)) != 2) stop("Please use exactly two unique groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few group labels.")
  group <- as.character(group)

  # Replace zeros unless alpha is provided
  if(missing(alpha)){ alpha <- as.numeric(NA)
  }else{ if(!is.na(alpha)) if(alpha == 0) alpha <- as.numeric(NA) }
  ct <- counts
  if(any(as.matrix(counts) == 0) & is.na(alpha)){
    message("Alert: Replacing 0s with next smallest value.")
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])
  }

  # Initialize @active, @weighted
  result <- new("propd")
  result@active <- "theta_d" # set theta_d active by default
  result@weights <- as.matrix(NA)
  result@weighted <- weighted
  result@dfz <- 0

  # Initialize @weights
  if(weighted){
    message("Alert: Calculating limma-based weights.")
    packageCheck("limma")
    design <- matrix(0, nrow = nrow(ct), ncol = 2)
    design[group == unique(group)[1], 1] <- 1
    design[group == unique(group)[2], 2] <- 1
    v <- limma::voom(t(counts), design = design)
    result@weights <- t(v$weights)
  }

  # Initialize @counts, @group, @alpha
  result@counts <- as.data.frame(ct)
  result@group <- as.character(group)
  result@alpha <- as.numeric(alpha)

  # Initialize @permutes
  result@permutes <- data.frame()
  if(p > 0){
    message("Alert: Fixing permutations to active random seed.")
    permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
    for(col in 1:ncol(permutes)) permutes[, col] <- sample(1:nrow(ct))
    result@permutes <- permutes
  }

  # Initialize @results
  result@results <-
    calculateTheta(result@counts, result@group, result@alpha,
                   weighted = result@weighted,
                   weights = result@weights)

  # Initialize @results -- Tally frequency of 0 counts
  # [note, use `counts` because `ct`` is already zero-replaced]
  if(any(as.matrix(counts) == 0)){
    message("Alert: Tabulating the presence of 0 counts.")
    result@results$Zeros <- ctzRcpp(as.matrix(counts)) # count 0s
  }

  # Initialize @results -- Round data to 14 digits
  if(any(result@results$theta > 1)){
    message("Alert: Theta rounded to 14 decimal digits.")
    result@results$theta <- round(result@results$theta, 14)
  }

  message("Alert: Use 'setActive' to select a theta type.")
  message("Alert: Use 'updateCutoffs' to calculate FDR.")
  message("Alert: Use 'updateF' to calculate F-stat.")

  return(result)
}

#' @rdname propd
#' @section Functions:
#' \code{setActive:}
#'  Build analyses and figures using a specific theta type. For
#'  example, set \code{what = "theta_d"} to analyze disjointed
#'  proportionality and \code{what = "theta_e"} to analyze
#'  emergent proportionality.
#' @export
setActive <- function(propd, what = "theta_d"){

  if(class(propd) != "propd") stop("Please provide a 'propd' object.")
  if(!any(what == colnames(propd@results)) & what != propd@active){
    stop("Provided theta type not recognized.")
  }

  # Rename old active theta type
  i <- which(colnames(propd@results) == "theta")
  colnames(propd@results)[i] <- propd@active

  # Name new active theta type
  i <- which(colnames(propd@results) == what)
  colnames(propd@results)[i] <- "theta"

  propd@active <- what
  message("Alert: Update FDR or F-stat manually.")

  return(propd)
}

#' @rdname propd
#' @section Functions:
#' \code{setDisjointed:}
#'  A wrapper for \code{setActive(propd, what = "theta_d")}.
#' @export
setDisjointed <- function(propd){
  setActive(propd, what = "theta_d")
}

#' @rdname propd
#' @section Functions:
#' \code{setEmergent:}
#'  A wrapper for \code{setActive(propd, what = "theta_e")}.
#' @export
setEmergent <- function(propd){

  if(table(propd@group)[1] != table(propd@group)[2]){

    warning("Emergent proportionality not yet validated for unequal group sizes.")
  }

  setActive(propd, what = "theta_e")
}
