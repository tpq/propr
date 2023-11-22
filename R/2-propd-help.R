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
#' @param object A \code{propd} object.
#' @param propd A \code{propd} object.
#'
#' @name propd
#' @importFrom methods show new
#' @importFrom Rcpp sourceCpp
NULL
