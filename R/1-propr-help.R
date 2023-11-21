#' The propr Package
#'
#' @description
#' Welcome to the \code{propr} package!
#'
#' To learn more about calculating proportionality, see
#'  Details.
#'
#' To learn more about differential proportionality, see
#'  \code{\link{propd}}.
#'
#' To learn more about compositional data analysis, see
#'  \code{citation("propr")}.
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
#' @param object A \code{propr} object.
#'
#' @name propr
#' @useDynLib propr, .registration = TRUE
#' @importFrom methods show new
#' @importFrom Rcpp sourceCpp
NULL
