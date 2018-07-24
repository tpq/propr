#' Build Index from ivar Argument
#'
#' This function builds an index from the \code{ivar} argument. Used by
#'  the \code{propr} initialize method and \code{updateF}.
#'
#' @inheritParams all
#' @return A numeric vector of indices for the reference set.
#'
#' @export
ivar2index <- function(counts, ivar){

  if(missing(ivar)) ivar <- 0
  if(!is.vector(ivar)) stop("Provide 'ivar' as vector.")
  `%is%` <- function(a, b) identical(a, b)

  if(ivar %is% 0 | ivar %is% NA | ivar %is% NULL | ivar %is% "all" | ivar %is% "clr"){

    use <- 1:ncol(counts) # use all features for geometric mean

  }else if(ivar %is% "iqlr"){

    if(any(counts == 0)){
      message("Alert: Replacing 0s in \"count matrix\" with 1 to calculate IQR.")
      counts[counts == 0] <- 1
    }

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

  return(use)
}

#' @rdname propr
#' @section Methods (by generic):
#' \code{subset:} Method to subset \code{propr} object.
#' @export
setMethod("subset", signature(x = "propr"),
          function(x, subset, select){

            if(missing(subset)) subset <- 1:nrow(x@counts)
            if(missing(select)) select <- 1:ncol(x@counts)

            if(is.character(select)){

              select <- match(select, colnames(x@counts))
            }

            x@counts <- x@counts[subset, select, drop = FALSE]
            x@logratio <- x@logratio[subset, select, drop = FALSE]
            x@matrix <- x@matrix[select, select, drop = FALSE]

            if(length(x@pairs) > 0){

              message("Alert: User must repopulate @pairs slot after `subset`.")
              x@pairs <- vector("numeric")
            }

            message("Alert: Using 'subset' is not compatible the @results table.")
            x@results <- data.frame()

            message("Alert: Using 'subset' disables permutation testing.")
            x@permutes <- list(NULL)

            return(x)
          }
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{[:} Method to subset \code{propr} object.
#' @aliases [,propr-method
#' @docType methods
#' @export
setMethod('[', signature(x = "propr", i = "ANY", j = "ANY"),
          function(x, i = "all", j, tiny = FALSE){

            if(i == "all"){

              x@pairs <- indexPairs(x@matrix, "all")
              return(x)
            }

            if(!i %in% c("==", "=", ">", ">=", "<", "<=", "!=")){

              stop("Operator not recognized. Index using e.g., `prop[\">\", .95]`.")
            }

            if(missing(j) | !is.numeric(j) | length(j) != 1){

              stop("Reference not found. Index using e.g., `prop[\">\", .95]`.")
            }

            newPairs <- indexPairs(x@matrix, i, j)
            if(length(newPairs) == 0) message("Alert: Method failed to index any pairs.")
            if(length(x@pairs) > 0) message("Alert: Appending prior index.")
            x@pairs <- sort(union(x@pairs, newPairs))

            if(tiny){

              return(simplify(x))

            }else{

              return(x)
            }
          }
)

#' @rdname propr
#' @section Functions:
#' \code{simplify:}
#'  This convenience function takes an indexed \code{propr} object
#'  and subsets the object based on that index. Then, it populates the
#'  \code{@@pairs} slot of the new object with an updated version
#'  of the original index. You can call \code{simplify} from within the
#'  \code{[} method using the argument \code{tiny}.
#' @export
simplify <- function(object){

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  # Subset propr object based on index
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  selection <- sort(union(coords[[1]], coords[[2]]))
  object@pairs <- vector("numeric")
  new <- subset(object, select = selection)

  # Repopulate the pairs slot
  for(i in 1:length(coords[[1]])){

    coords[[1]][i] <- which(selection == coords[[1]][i])
    coords[[2]][i] <- which(selection == coords[[2]][i])
  }

  new@pairs <- (coords[[2]] - 1) * nrow(new@matrix) + (coords[[1]] - 1) + 1

  return(new)
}

#' @rdname propr
#' @section Functions:
#' \code{updateCutoffs:}
#'  Use the \code{propr} object to permute proportionality
#'  across a number of cutoffs. Since the permutations get saved
#'  when the object is created, calling \code{updateCutoffs}
#'  will use the same random seed each time.
#' @export
updateCutoffs.propr <- function(object, cutoff = seq(.05, .95, .3)){

  if(identical(object@permutes, list(NULL))) stop("Permutation testing is disabled.")

  # Let NA cutoff skip function
  if(identical(cutoff, NA)) return(object)

  # Set up FDR cutoff table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
  FDR$cutoff <- cutoff
  p <- length(object@permutes)

  # Calculate propr for each permutation -- NOTE: `select` and `subset` disable permutation testing
  for(k in 1:p){

    numTicks <- progress(k, p, numTicks)

    # Calculate propr exactly based on @metric, @ivar, and @alpha
    ct.k <- object@permutes[[k]]
    pr.k <- suppressMessages(
      propr(ct.k, object@metric, ivar = object@ivar, alpha = object@alpha))
    pkt <- pr.k@results$propr

    # Find number of permuted theta less than cutoff
    for(cut in 1:nrow(FDR)){ # randcounts as cumsum

      # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
      if(object@metric == "rho" | object@metric == "cor"){
        FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + sum(pkt > FDR[cut, "cutoff"])
      }else{ # phi & phs
        FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + sum(pkt < FDR[cut, "cutoff"])
      }
    }
  }

  # Calculate FDR based on real and permuted tallys
  FDR$randcounts <- FDR$randcounts / p # randcounts as mean
  for(cut in 1:nrow(FDR)){

    # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
    if(object@metric == "rho" | object@metric == "cor"){
      FDR[cut, "truecounts"] <- sum(object@results$propr > FDR[cut, "cutoff"])
    }else{ # phi & phs
      FDR[cut, "truecounts"] <- sum(object@results$propr < FDR[cut, "cutoff"])
    }

    FDR[cut, "FDR"] <- FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
  }

  # Initialize @fdr
  object@fdr <- FDR

  return(object)
}
