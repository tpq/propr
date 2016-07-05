#' Calculate proportionality metric rho.
#'
#' \code{perb} returns a \code{propr} object containing measures of proportionality.
#'
#' Calculates proportionality metric rho described in Lovell 2015 and expounded
#'  in Erb 2016. Uses centered log-ratio transformation of data by default,
#'  but will use additive log-ratio transformation of data if non-zero
#'  \code{ivar} provided.
#'
#' Let d represent any number of features measured across multiple biological replicates n
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with n rows and d columns into a proportionality matrix of d rows and d
#' 	columns containing rho measurements for each feature pair. One can think of the resultant
#' 	matrix as equivalent to a correlation matrix.
#'
#' @param ivar A numeric scalar. Specificies feature to use as reference for additive log-ratio transformation.
#' @inheritParams phit
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{phit}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' prop <- perb(counts, ivar = 0, iter = 0)
#' @importFrom methods new
#' @importFrom stats ecdf p.adjust
#' @export
perb <- function(counts, ivar = 0, iter = 0, iterSize = ncol(counts) - (ivar > 0), iterHow = 1, onlyDistr = FALSE){

  if(!onlyDistr){

    cat("Calculating all rho for \"count matrix\"...\n")
    prop <- new("propr")
    prop@counts <- as.data.frame(counts)
    if(ivar != 0){ prop@logratio <- as.data.frame(proprALR(prop@counts, ivar))
    }else{ prop@logratio <- as.data.frame(proprCLR(prop@counts))}
    prop@matrix <- proprPerb(prop@counts, ivar)
    prop@pairs <- proprPairs(prop@matrix)

    prop@pairs <- prop@pairs[rev(order(abs(prop@pairs$prop))), ]
    rownames(prop@pairs) <- 1:nrow(prop@pairs)
    if(iter == 0) return(prop)

  }else{

    if(iter == 0){

      stop("This function cannot return a fit distribution if iter = 0!")
    }
  }

  distr <- vector("numeric", iter * iterSize * (iterSize - 1) / 2)
  for(i in 1:iter){

    cat(paste0("Calculating simulated rho for iter ", i, "...\n"))

    # Handle alr properly
    if(ivar != 0){

      if(iterHow == 1){

        index.i <- sample(1:ncol(counts[, -ivar]), iterSize)
        counts.i <- apply(counts[, -ivar][, index.i], 2, sample)

      }else if(iterHow == 2){


        counts.i <- t(apply(counts[, -ivar], 1, sample, iterSize))

      }else{

        stop("Uh oh! Provided 'iterHow' not recognized! Select '1' for features and '2' for subject.")
      }

      fixed <- counts[, ivar]
      counts.i <- cbind(counts.i, fixed)
      ivar.i <- ncol(counts.i)

    }else{

      if(iterHow == 1){

        index.i <- sample(1:ncol(counts), iterSize)
        counts.i <- apply(counts[, index.i], 2, sample)

      }else if(iterHow == 2){

        counts.i <- t(apply(counts, 1, sample, iterSize))

      }else{

        stop("Uh oh! Provided 'iterHow' not recognized! Select '1' for features and '2' for subject.")
      }

      ivar.i <- 0
    }

    prop.i <- proprPerb(counts.i, ivar.i)
    begin <- (i - 1) * iterSize * (iterSize - 1) / 2 + 1
    end <- (i - 1) * iterSize * (iterSize - 1) / 2 + iterSize * (iterSize - 1) / 2
    distr[begin:end] <- proprTri(prop.i)
    rm(prop.i)
  }

  cat("Fitting rho to distribution...\n")
  fit <- ecdf(1 - abs(distr))
  rm(distr)

  if(!onlyDistr){

    cat("Using fit to convert rho into pval...\n")
    pval <- fit(1 - abs(prop@pairs$prop))
    # pval[pval >= .5] <- 1 - pval[pval >= .5] # make 2-tails
    # pval <- pval * 2 # scale 0 to 1

    cat("Correcting for multiple testing...\n")
    fdr <- p.adjust(pval, method = "BH")

    cat("Building results...\n")
    prop@pairs <- data.frame(prop@pairs, "pval" = pval, "fdr" = fdr, stringsAsFactors = FALSE)
    # prop@pairs <- prop@pairs[order(prop@pairs$pval), ]

    return(prop)

  }else{

    return(fit)
  }
}
