#' Calculate Proportionality
#'
#' @description
#' Let D represent any number of features measured across N biological replicates
#' 	exposed to a binary or continuous event E. For example, E could indicate case-control
#' 	status, treatment status, treatment dose, or time. This function converts a "count matrix"
#' 	with N rows and D columns into a proportionality matrix of D rows and D columns.
#'
#' For phi, the result of \code{phit}, one can think of the resultant matrix as
#' 	analogous to a distance matrix, except that it has no symmetry unless forced.
#' 	For rho, the result of \code{perb}, one can think of the resultant matrix as
#' 	analogous to a correlation matrix.
#' 	For phs, the result of \code{phis}, one can think of the resultant matrix as
#' 	either a naturally symmetric variant of phi or a monotonic variant of rho.
#'
#' These methods all use the centered log-ratio transformation by default,
#'  but will use an additive log-ratio transformation instead if a scalar
#'  \code{ivar} is provided. When using an additive log-ratio transformation,
#'  this function will return \code{rho = 0} for the column and row in the
#'  \code{@@matrix} slot that would contain the reference feature.
#'  Setting \code{ivar} to a numeric or character vector will transform
#'  data using the geometric mean of only the indexed features.
#'  Alternatively, setting \code{ivar} to "iqlr" will transform data using
#'  the geometric mean of only the features with variances that fall in
#'  the inter-quartile range of all per-feature variances. We base this
#'  "iqlr" transformation on the \code{ALDEx2} package.
#'
#' Log-ratio transformation, by its nature, fails if the input data
#'  contain any zero values. To avoid an error in this case, these
#'  methods automatically replace all zero values with 1. However,
#'  the topic of zero replacement is controversial. Proceed carefully
#'  when analyzing data that contain any zero values.
#'
#' The \code{select} argument subsets the feature matrix
#'  after log-ratio transformation but before calculating
#'  proportionality. This reduces the run-time and RAM
#'  overhead without impacting the final result. Removing
#'  lowly abundant features prior to log-ratio transformation
#'  could otherwise change the proportionality measure.
#'
#' @param counts A data.frame or matrix. A "count matrix" with
#'  subjects as rows and features as columns.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry
#'  by reflecting the "lower left triangle".
#' @param ivar A numeric scalar. Specificies reference feature(s)
#'  for additive log-ratio transformation. The argument will also
#'  accept feature name(s) instead of the index position(s).
#'  Set to "iqlr" to use inter-quartile log-ratio transformation.
#'  Ignore to use centered log-ratio transformation.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Optional. Use this argument to subset the proportionality
#'  matrix before building without altering the final result.
#' @param .Object Missing. Ignore. Leftover from the generic
#'  method definition.
#'
#' @return Returns a \code{propr} object.
#'
#' @examples
#' library(propr)
#' data(mail)
#' phi <- phit(mail)
#' rho <- perb(mail)
#' phs <- phis(mail)
#' @name proportionality
NULL

#' @rdname proportionality
setMethod("initialize", "propr",
          function(.Object, counts, ivar, select){

            # Retain backwards-compatibility
            if(missing(counts)){
              return(.Object)
            }

            # Quality control check input
            if(any(is.na(counts))) stop("Remove NAs from 'counts' before proceeding.")
            if(class(counts) == "data.frame") counts <- as.matrix(counts)
            if(is.null(colnames(counts))) colnames(counts) <- as.character(1:ncol(counts))
            if(is.null(rownames(counts))) rownames(counts) <- as.character(1:nrow(counts))

            # Replace zeros if needed
            if(any(0 == counts)){
              message("Alert: Replacing 0s in \"count matrix\" with 1.")
              counts[counts == 0] <- 1
            }

            # Get feature set to use in g(x) calculation
            if(missing(ivar)) ivar <- 0
            if(!is.vector(ivar)) stop("Provide 'ivar' as vector.")
            `%is%` <- function(a, b) identical(a, b)
            if(ivar %is% 0 | ivar %is% NA | ivar %is% NULL | ivar %is% "all" | ivar %is% "clr"){

              use <- 1:ncol(counts) # use all features for geometric mean

            }else if(ivar %is% "iqlr"){

              counts.clr <- apply(log(counts), 1, function(x){ x - mean(x) })
              counts.var <- apply(counts.clr, 1, var)
              quart <- stats::quantile(counts.var) # use features with unextreme variance
              use <- which(counts.var < quart[4] & counts.var > quart[2])

            }else{

              if(is.character(ivar)){
                if(!all(ivar %in% colnames(counts))) stop("Some 'ivar' not in data.")
                use <- which(colnames(counts) %in% ivar) # use features given by name
              }else{
                use <- sort(ivar) # use features given by number
              }
            }

            # Use g(x) to log-ratio transform data
            logX <- log(counts)
            logSet <- logX[, use, drop = FALSE]
            lr <- sweep(logX, 1, rowMeans(logSet), "-")

            # Subset data after transformation
            if(!missing(select)){

              # Make select boolean (it's OK if it's integer)
              if(!is.vector(select)) stop("Provide 'select' as vector.")
              if(is.character(select)) select <- match(select, colnames(counts))
              if(any(is.na(select))) stop("Some 'select' not in data.")
              counts <- counts[, select]
              lr <- lr[, select]
            }

            .Object@counts <- as.matrix(counts)
            .Object@logratio <- lr
            .Object@pairs <- vector("numeric")

            return(.Object)
          }
)

#' @rdname proportionality
#' @export
phit <- function(counts, ivar = 0, select, symmetrize = TRUE){

  prop <- new("propr", counts = counts, ivar = ivar, select = select)
  prop@matrix <- lr2phi(prop@logratio)
  if(symmetrize) symRcpp(prop@matrix)
  return(prop)
}

#' @rdname proportionality
#' @export
perb <- function(counts, ivar = 0, select){

  prop <- new("propr", counts = counts, ivar = ivar, select = select)
  prop@matrix <- lr2rho(prop@logratio)
  return(prop)
}

#' @rdname proportionality
#' @export
phis <- function(counts, ivar = 0, select){

  prop <- new("propr", counts = counts, ivar = ivar, select = select)
  prop@matrix <- lr2phs(prop@logratio)
  return(prop)
}
