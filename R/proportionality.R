#' Calculate Proportionality
#'
#' @description
#' Let D represent any number of features measured across N biological replicates
#' 	exposed to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a "count matrix"
#' 	with N rows and D columns into a proportionality matrix of D rows and D columns.
#'
#' For phi, the result of \code{phit}, one can think of the resultant matrix as
#' 	analogous to a distance matrix, except that it has no symmetry unless forced.
#' 	For phs, the result of \code{phis}, one can think of the resultant matrix as
#' 	either a naturally symmetric variant of phi or a monotonic variant of rho.
#' 	For rho, the result of \code{perb}, one can think of the resultant matrix as
#' 	analogous to a correlation matrix.
#'
#' \code{perb} and \code{phis} use a centered log-ratio transformation by default,
#'  but will use an additive log-ratio transformation instead if a non-zero
#'  \code{ivar} is provided. When using an additive log-ratio transformation,
#'  this function will return \code{rho = 0} for the column and row in the
#'  \code{@@matrix} slot that would contain the reference feature.
#'
#' Log-ratio transformation, by its nature, fails if the input data contain
#'  any zero values. To avoid an error in this case, these functions automatically
#'  replace all zero values with 1. Note, however, that the topic of
#'  zero replacement is controversial. Proceed carefully when analyzing data
#'  that contain zero values.
#'
#' @param counts A data.frame or matrix. A "count matrix" with
#'  subjects as rows and features as columns.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry
#'  by reflecting the "lower left triangle".
#' @param ivar A numeric scalar. Specificies a reference feature
#'  for additive log-ratio transformation. The argument will also
#'  accept a feature name instead of the index position.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Optional. Use this argument to subset the proportionality
#'  matrix without altering the final result.
#'
#' @return Returns a \code{propr} object.
#'
#' @examples
#' library(propr)
#' data(mail)
#' phi <- phit(mail)
#' phs <- phis(mail)
#' rho <- perb(mail)
#' @name proportionality
NULL

#' @rdname proportionality
#' @export
phit <- function(counts, symmetrize = TRUE){

  if(any(is.na(counts))) stop("Uh oh! Remove NAs before proceeding.")
  prop <- new("propr")
  prop@counts <- as.matrix(counts)

  if(any(0 == prop@counts)){

    message("Alert: Replacing 0s in \"count matrix\" with 1.")
    prop@counts[prop@counts == 0] <- 1
  }

  prop@logratio <- clrRcpp(prop@counts[]) # [] forces copy
  prop@matrix <- phiRcpp(prop@counts[], symmetrize) # [] forces copy
  prop@pairs <- vector("numeric")

  return(prop)
}

#' @rdname proportionality
#' @export
perb <- function(counts, ivar = 0, select){

  if(any(is.na(counts))) stop("Uh oh! Remove NAs before proceeding.")
  prop <- new("propr")
  prop@counts <- as.matrix(counts)

  if(any(0 == prop@counts)){

    message("Alert: Replacing 0s in \"count matrix\" with 1.")
    prop@counts[prop@counts == 0] <- 1
  }

  if(ivar != 0){

    if(is.character(ivar)){

      # Find i-th index of ivar name
      index <- ivar == colnames(prop@counts)
      if(!any(index)) stop("Uh oh! Provided ivar reference not found in data.")
      ivar <- which(index)
    }

    prop@logratio <- alrRcpp(prop@counts[], ivar)

  }else{

    prop@logratio <- clrRcpp(prop@counts[])
  }

  if(!missing(select)){

    # Make select boolean (it's OK if it's integer)
    if(is.character(select)) select <- match(select, colnames(prop@counts))
    if(any(is.na(select))) stop("Uh oh! Provided select reference not found in data.")

    # Map ivar to new subset (else assign it 0)
    mapping <- (1:ncol(prop@counts))[select]
    if(any(ivar == mapping)){ ivar <- which(ivar == mapping)
    }else{ ivar <- 0}

    # Now OK to drop the unchanged reference
    prop@counts <- prop@counts[, select]
    prop@logratio <- prop@logratio[, select]
  }

  prop@matrix <- rhoRcpp(prop@counts[], prop@logratio[], ivar)
  prop@pairs <- vector("numeric")

  return(prop)
}

#' @rdname proportionality
#' @export
phis <- function(counts, ivar = 0, select){

  prop <- perb(counts, ivar, select)
  rhoToPhs(prop@matrix)
  return(prop)
}
