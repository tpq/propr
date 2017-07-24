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
#' The \code{propd} method calculates theta. However, VLR fails in
#'  the setting of zero counts. The \code{propd} method
#'  will use a Box-Cox transformation to approximate VLR based on
#'  the parameter \eqn{\alpha} if provided. We refer the user to
#'  the vignette for more details.
#'
#' Note that Group 1 always refers to the first element of the
#'  \code{group} vector argument supplied to \code{propd}.
#'
#' @details
#' \code{setActive:}
#'  Build analyses and figures using a specific theta type. For
#'  example, set \code{what = "theta_d"} to analyze disjointed
#'  proportionality and \code{what = "theta_e"} to analyze
#'  emergent proportionality. These provide the same results as
#'  \code{setDisjointed} and \code{setEmergent}, respectively.
#'
#' \code{updateCutoffs:}
#'  Use the \code{propd} object to permute theta across a
#'  number of theta cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#'
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
#'
#' \code{shale:}
#'  Builds a table of within-group and total log-ratio
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
#' \code{geyser:}
#'  Plots indexed pairs based on the within-group
#'  log-ratio variance (VLR) for each group. Pairs near the
#'  origin show a highly proportional relationship in
#'  both groups. Each line away from the \code{y = x} line
#'  indicates a doubling of VLR compared to the other group.
#'  All pairs colored based on PAL
#'  (see: \code{\link{pals}}).
#'  See \code{gemini}.
#'
#' \code{bowtie:}
#'  Plots indexed pairs based on the log-ratio means
#'  (LRM), relative to its PAL, for each group. Pairs near
#'  the origin show comparable LRM, relative to its PAL, in
#'  both groups. Each line away from the \code{y = x} line
#'  indicates a doubling of LRM compared to the other group.
#'  All pairs colored based on PAL
#'  (see: \code{\link{pals}}).
#'  See \code{gemini}.
#'
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
#'
#' \code{slice:}
#'  Plots log-ratio counts relative to a \code{reference}
#'  node for all pairs that include the reference itself.
#'  This makes a useful adjunct function to visualize how
#'  features vary across samples relative to a PAL.
#'
#' \code{decomposed:}
#'  Plots the decomposition of log-ratio variance into
#'  (weighted) group variances and between-group variance.
#'  Useful for visualizing how a theta type selects pairs.
#'
#' \code{propd2propr:}
#'  Transforms a \code{propd} object into a \code{propr} object
#'  where the \code{@@matrix} slot contains \eqn{1 - \theta}.
#'  Allows the user to interrogate theta using any
#'  function described in \code{\link{visualize}}.
#'
#' @slot counts A data.frame. Stores the original "count matrix" input.
#' @slot group A character vector. Stores the original group labels.
#' @slot alpha A double. Stores the alpha value used for transformation.
#' @slot weighted A logical. Stores whether the theta is weighted.
#'  Used by \code{updateCutoffs}.
#' @slot active A character. Stores the name of the active theta type.
#'  Used by \code{updateCutoffs}.
#' @slot theta A data.frame. Stores the pairwise theta measurements.
#' @slot permutes A data.frame. Stores the shuffled group labels,
#'  used to reproduce permutations of theta.
#' @slot fdr A data.frame. Stores the FDR cutoffs for theta.
#'
#' @param object,x,propd A \code{propd} object.
#' @param y Missing. Ignore. Leftover from the generic method
#'  definition.
#' @param counts A data.frame or matrix. A "count matrix" with
#'  subjects as rows and features as columns.
#' @param ivar A numeric scalar. Specifies reference feature(s)
#'  for additive log-ratio transformation. The argument will also
#'  accept feature name(s) instead of the index position(s).
#'  Set to "iqlr" to use inter-quartile log-ratio transformation.
#'  Ignore to use centered log-ratio transformation.
#' @param what A character string. The theta type to set active.
#' @param group A character vector. Group or sub-group memberships,
#'  ordered according to the row names in \code{counts}.
#' @param alpha A double. See vignette for details. Leave missing
#'  to skip Box-Cox transformation.
#' @param weighted A boolean. Toggles whether to calculate
#'  theta using \code{limma::voom} weights.
#' @param p An integer. The number of permutation cycles.
#' @param cutoff For \code{propd}, a numeric vector.
#'  this argument provides the FDR cutoffs to test for theta.
#'  For graph functions, a numeric scalar. This argument
#'  indicates the maximum theta to include in the figure.
#'  For graph functions, a large integer will instead
#'  retrieve the top N pairs as ranked by theta.
#' @param k An integer. The maximum number of PALs to index
#'  when calculating \code{\link{pals}} in the network.
#' @param reference A character string. A feature to use as the
#'  denominator reference when comparing log-ratio abundances.
#' @param propr An indexed \code{propr} object. Use to add
#'  proportional edges (colored green) to a \code{propd} network.
#' @param clean A boolean. Toggles whether to remove pairs
#'  with "Bridged" or "Missing" PALs. Used by \code{geyser},
#'  \code{bowtie}, and \code{gemini}.
#' @param plotSkip A boolean. Toggles whether to build
#'  the network graph without plotting it.
#'  Used by \code{\link{pals}}.
#' @inheritParams visualize
#'
#' @name propd
#' @importFrom Rcpp sourceCpp
#' @importFrom methods show new
#' @importFrom utils head
#' @importFrom stats var
NULL

#' @rdname propd
#' @export
setClass("propd",
         slots = c(
           counts = "data.frame",
           group = "character",
           alpha = "numeric",
           weighted = "logical",
           active = "character",
           theta = "data.frame",
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

            cat("@counts summary:",
                nrow(object@counts), "subjects by", ncol(object@counts), "features\n")

            cat("@group summary:", length(unique(object@group)), "unique groups (",
                paste(table(object@group), collapse = " x "), ")\n")

            cat("@theta summary:",
                nrow(object@theta), "feature pairs (", object@active, ")\n")

            cat("@fdr summary:",
                ncol(object@permutes), "iterations\n")

            print(object@fdr)

            cat("See ?propd for object methods\n")
          }
)

#' @rdname propd
#' @export
propd <- function(counts, group, alpha, p = 100, cutoff = NA,
                  weighted = FALSE){

  # Clean "count matrix"
  if(any(is.na(counts))) stop("Remove NAs from 'counts' before proceeding.")
  if(class(counts) == "data.frame") counts <- as.matrix(counts)
  if(is.null(colnames(counts))) colnames(counts) <- as.character(1:ncol(counts))
  if(is.null(rownames(counts))) rownames(counts) <- as.character(1:nrow(counts))
  if(missing(alpha)) alpha <- NA
  ct <- counts

  # Replace zeros unless alpha is provided
  if(any(as.matrix(counts) == 0) & is.na(alpha)){
    message("Alert: Replacing 0s in \"count matrix\" with 1.")
    ct[ct == 0] <- 1
  }

  # Prepare theta results object
  result <- new("propd")
  result@theta <- calculateTheta(ct, group, alpha, weighted = weighted)
  result@weighted <- weighted
  result@active <- "theta_d" # set theta_d active by default

  # Tally frequency of 0 counts
  if(any(as.matrix(counts) == 0)){
    message("Alert: Tabulating the presence of 0 counts.")
    result@theta$Zeros <- ctzRcpp(as.matrix(counts)) # count 0s before replacement
  }

  # Save important intermediate values
  result@counts <- as.data.frame(ct)
  result@group <- as.character(group)
  if(!missing(alpha)){ result@alpha <- as.numeric(alpha)
  }else{ result@alpha <- as.numeric(NA) }

  # Pre-compute group samples for permutation
  permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
  for(col in 1:ncol(permutes)) permutes[, col] <- sample(1:nrow(ct))
  result@permutes <- permutes

  # Round data to 14 decimal places
  if(any(result@theta$theta > 1)){
    message("Alert: Theta rounded to 14 decimal digits.")
    result@theta$theta <- round(result@theta$theta, 14)
    result@theta$theta_e <- round(result@theta$theta_e, 14)
    result@theta$theta_f <- round(result@theta$theta_f, 14)
  }

  # Calculate FDR for cutoffs
  result <- updateCutoffs(result, cutoff)

  return(result)
}

#' @rdname propd
#' @export
propd2propr <- function(object, ivar){

  prop <- new("propr", object@counts, ivar)
  prop@matrix <- 1 - half2mat(object@theta$theta)
  return(prop)
}

#' @rdname propd
#' @export
setActive <- function(propd, what = "theta_d"){

  if(class(propd) != "propd") stop("Please provide a 'propd' object.")
  if(!any(what == colnames(propd@theta)) & what != propd@active){
    stop("Provided theta type not recognized.")
  }

  # Rename old active theta type
  i <- which(colnames(propd@theta) == "theta")
  colnames(propd@theta)[i] <- propd@active

  # Name new active theta type
  i <- which(colnames(propd@theta) == what)
  colnames(propd@theta)[i] <- "theta"

  propd@active <- what
  message("Use 'updateCutoffs' to refresh FDR.")
  return(propd)
}

#' @rdname propd
#' @export
setDisjointed <- function(propd){

  setActive(propd, what = "theta_d")
}

#' @rdname propd
#' @export
setEmergent <- function(propd){

  setActive(propd, what = "theta_e")
}
