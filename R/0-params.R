#' A list of most parameters.
#'
#' @description
#' A list of most parameters.
#'
#' @param counts A data.frame or matrix. A "count matrix" with
#'  subjects as rows and features as columns. Note that this matrix
#'  does not necessarily have to contain counts.
#' @param metric A character string. The proportionality metric
#'  to calculate. Choose from "rho", "phi", or "phs".
#' @param ivar A numeric scalar. Specifies reference feature(s)
#'  for additive log-ratio transformation. The argument will also
#'  accept feature name(s) instead of the index position(s).
#'  Set to "iqlr" to use inter-quartile log-ratio transformation.
#'  Ignore to use centered log-ratio transformation.
#' @param select Optional. Use this to subset the final
#'  proportionality matrix without altering the result.
#'  Use this argument to rearrange feature order.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry
#'  by reflecting the "lower left triangle".
#' @param alpha A double. See vignette for details. Leave missing
#'  to skip Box-Cox transformation.
#' @param p An integer. The number of permutation cycles.
#' @param ... Arguments passed to the wrapped method.
#'
#' @param object,x,rho,propr,propd A \code{propr} or \code{propd} object.
#' @param subset Subsets via \code{object@counts[subset, ]}.
#'  Use this argument to rearrange subject order.
#'  For backwards compatibility.
#' @param i Operation used for the subset indexing. Select from
#'  "==", "=", ">", ">=", "<", "<=", "!=", or "all".
#'  For backwards compatibility.
#' @param j Provide a numeric value to which to compare the
#'  proportionality measures in the \code{@@matrix} slot.
#'  For backwards compatibility.
#' @param tiny A logical scalar. Toggles whether to pass the indexed
#'  result through \code{\link{simplify}}.
#'  For backwards compatibility.
#'
#' @param group A character vector. Group or sub-group memberships,
#'  ordered according to the row names in \code{counts}.
#' @param weighted A boolean. Toggles whether to calculate
#'  theta using \code{limma::voom} weights.
#' @param cutoff For \code{updateCutoffs}, a numeric vector.
#'  this argument provides the FDR cutoffs to test.
#'  For graph functions, a numeric scalar. This argument
#'  indicates the maximum theta to include in the figure.
#'  For graph functions, a large integer will instead
#'  retrieve the top N pairs as ranked by theta.
#' @param what A character string. The theta type to set active.
#' @param moderated For \code{updateF}, a boolean. Toggles
#'  whether to calculate a moderated F-statistic.
#'
#' @param y Missing. Ignore. Leftover from the generic
#'  method definition.
#' @param prompt A logical scalar. Set to \code{FALSE} to disable
#'  the courtesy prompt when working with big data.
#' @param plotly A logical scalar. Set to \code{TRUE} to produce
#'  a dynamic plot using the \code{plotly} package.
#' @param k An integer. For \code{propr} methods, the number
#'  of co-clusters (where all pairs receive a specified color
#'  if and only if both members belong to same the cluster).
#'  For \code{propd} methods, the maximum number of PALs to index
#'  when calculating \code{\link{pals}} in the network.
#' @param col1,col2 A character vector. Specifies which nodes
#'  to color \code{red} or \code{blue}, respectively.
#' @param d3 A boolean. Use \code{rgl} to plot 3D network.
#'
#' @param include This argument indicates which features by
#'  name should belong to a pair for that pair to get included
#'  in the results. Subset performed by
#'  \code{Partner \%in\% subset | Pair \%in\% subset}.
#' @param or A boolean. If \code{FALSE}, \code{include} subsets
#'  by \code{Partner \%in\% subset & Pair \%in\% subset}.
#' @param clean A boolean. Toggles whether to remove pairs
#'  with "Bridged" or "Missing" PALs. Used by \code{geyser},
#'  \code{bowtie}, and \code{gemini}.
#' @param plotSkip A boolean. Toggles whether to build
#'  the network graph without plotting it.
#'  Used by \code{\link{pals}}.
#'
#' @name all
NULL
