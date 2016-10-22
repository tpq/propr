#' Calculate proportionality metric rho.
#'
#' \code{perb} returns a \code{propr} object containing a measure of proportionality.
#'
#' Let d represent any number of features measured across n biological replicates
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with n rows and d columns into a proportionality matrix of d rows and d
#' 	columns containing rho measurements for each feature pair. One can think of the resultant
#' 	matrix as equivalent to a correlation matrix.
#'
#' This function uses a centered log-ratio transformation of the data by default,
#'  but will use an additive log-ratio transformation of the data if a non-zero
#'  \code{ivar} is provided. When using an additive log-ratio transformation,
#'  this function will return \code{rho = 0} for each pair containing the
#'  reference feature.
#'
#' @param ivar A numeric scalar. Specificies a reference feature
#'  for additive log-ratio transformation. Argument will now accept
#'  the feature name instead of index.
#' @param select Subsets via \code{object@counts[, select]}.
#'  Optional. Use this argument to subset the proportionality
#'  matrix without impacting the value of \code{rho}.
#' @inheritParams phit
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{phit}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' rho <- perb(counts, ivar = 0)
#' @importFrom methods new
#' @export
perb <- function(counts, ivar = 0, select){

  cat("Calculating rho from \"count matrix\".\n")
  prop <- new("propr")
  prop@counts <- as.matrix(counts)

  if(any(0 == prop@counts)){

    message("Replacing zeroes in \"count matrix\" with next smallest value.")
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
