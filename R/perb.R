#' Calculate phi
#'
#' \code{phit} returns a \code{propr} object containing a measure of proportionality, phi.
#'
#' Let D represent any number of features measured across N biological replicates
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with N rows and D columns into a proportionality matrix of D rows and D
#' 	columns containing a value of phi for each feature pair. One can think of the resultant
#' 	matrix as analogous to a distance matrix, except that it has no symmetry unless forced.
#'
#' @param counts A data.frame or matrix. A "count matrix" with subjects as rows and features as columns.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry by reflecting the "lower left triangle".
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{perb}}
#'
#' @examples
#' library(propr)
#' data(mail)
#' phi <- phit(mail, symmetrize = TRUE)
#' @importFrom methods new
#' @export
phit <- function(counts, symmetrize = TRUE){

  prop <- new("propr")
  prop@counts <- as.matrix(counts)

  if(any(0 == prop@counts)){

    message("Alert: Replacing 0s in \"count matrix\" with next smallest value.")
    prop@counts[prop@counts == 0] <- unique(sort(prop@counts))[2]
  }

  prop@logratio <- clrRcpp(prop@counts[]) # [] forces copy
  prop@matrix <- phiRcpp(prop@counts[], symmetrize) # [] forces copy
  prop@pairs <- vector("numeric")

  return(prop)
}

#' Calculate rho
#'
#' \code{perb} returns a \code{propr} object containing a measure of proportionality, rho.
#'
#' Let D represent any number of features measured across N biological replicates
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with N rows and D columns into a proportionality matrix of D rows and D
#' 	columns containing a value of rho for each feature pair. One can think of the resultant
#' 	matrix as analogous to a correlation matrix.
#'
#' This function uses a centered log-ratio transformation of the data by default,
#'  but will use an additive log-ratio transformation instead if a non-zero
#'  \code{ivar} is provided. When using an additive log-ratio transformation,
#'  this function will return \code{rho = 0} for the column and row in the
#'  \code{@@matrix} slot that would contain the reference feature.
#'
#' @param ivar A numeric scalar. Specificies a reference feature
#'  for additive log-ratio transformation. The argument will also
#'  accept a feature name instead of the index position.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Optional. Use this argument to subset the proportionality
#'  matrix without altering the final value of \code{rho}.
#' @inheritParams phit
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{phit}}
#'
#' @examples
#' library(propr)
#' data(mail)
#' rho <- perb(mail, ivar = 0)
#' @importFrom methods new
#' @export
perb <- function(counts, ivar = 0, select){

  prop <- new("propr")
  prop@counts <- as.matrix(counts)

  if(any(0 == prop@counts)){

    message("Alert: Replacing 0s in \"count matrix\" with next smallest value.")
    prop@counts[prop@counts == 0] <- unique(sort(prop@counts))[2]
  }

  if(class(ivar) == "character"){

    index <- ivar == colnames(prop@counts)
    if(any(index)){ ivar <- which(index)
    }else{ stop("Uh oh! Provided ivar reference not found in data.")}
  }

  if(ivar != 0){ prop@logratio <- alrRcpp(prop@counts[], ivar) # [] forces copy
  }else{ prop@logratio <- clrRcpp(prop@counts[])} # [] forces copy

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
