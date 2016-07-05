#' Calculate proportionality metric phi.
#'
#' \code{phit} returns a \code{propr} object containing measures of proportionality.
#'
#' Let d represent any number of features measured across multiple biological replicates n
#' 	subjected to a binary or continuous event E. For example, E could represent case-control
#' 	status, treatment status, treatment dose, or time. This function converts a
#' 	"count matrix" with n rows and d columns into a proportionality matrix of d rows and d
#' 	columns containing phi measurements for each feature pair. One can think of the resultant
#' 	matrix as equivalent to a distance matrix, except that it has no symmetry by default.
#'
#' @param counts A data.frame or matrix. A "count matrix" with subjects as rows and features as columns.
#' @param symmetrize A logical. If \code{TRUE}, forces symmetry by duplicating the "lower left triangle".
#' @param iter A numeric scalar. Fits \code{iter*iterSize*(iterSize-1)/2} values to an empiric distribution. Skip with \code{iter = 0}.
#' @param iterSize A numeric scalar. Fits \code{iter*iterSize*(iterSize-1)/2} values to an empiric distribution.
#' @param iterHow A numeric scalar. Select \code{1} to randomize feature vectors or \code{2} to randomize subject vectors.
#' @param onlyDistr A logical. Provided for backend use. Evokes function to return only \code{ecdf} fit.
#' @return Returns a \code{propr} object.
#'
#' @seealso \code{\link{propr}}, \code{\link{propr-class}}, \code{\link{perb}}
#'
#' @examples
#' randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
#' counts <- matrix(randomNum, nrow = 25, ncol = 10)
#' prop <- phit(counts, symmetrize = TRUE, iter = 0)
#' @importFrom methods new
#' @importFrom stats ecdf p.adjust
#' @export
phit <- function(counts, symmetrize = TRUE, iter = 0, iterSize = ncol(counts), iterHow = 1, onlyDistr = FALSE){

  if(!onlyDistr){

    cat("Calculating all phi for \"count matrix\"...\n")
    prop <- new("propr")
    prop@counts <- as.data.frame(counts)
    prop@logratio <- as.data.frame(proprCLR(prop@counts))
    prop@matrix <- proprPhit(prop@counts, symmetrize)
    prop@pairs <- proprPairs(prop@matrix)
    if(iter == 0) return(prop)

  }else{

    if(iter == 0){

      stop("This function cannot return a fit distribution if iter = 0!")
    }
  }

  distr <- vector("numeric", iter * iterSize * (iterSize - 1) / 2)
  for(i in 1:iter){

    cat(paste0("Calculating simulated phi for iter ", i, "...\n"))

    if(iterHow == 1){

      index.i <- sample(1:ncol(counts), iterSize)
      counts.i <- apply(counts[, index.i], 2, sample)

    }else if(iterHow == 2){

      counts.i <- t(apply(counts, 1, sample, iterSize))

    }else{

      stop("Uh oh! Provided 'iterHow' not recognized! Select '1' for features and '2' for subject.")
    }

    prop.i <- proprPhit(counts.i)
    begin <- (i - 1) * iterSize * (iterSize - 1) / 2 + 1
    end <- (i - 1) * iterSize * (iterSize - 1) / 2 + iterSize * (iterSize - 1) / 2
    distr[begin:end] <- proprTri(prop.i)
    rm(prop.i)
  }

  cat("Fitting phi to distribution...\n")
  fit <- ecdf(distr)
  rm(distr)

  if(!onlyDistr){

    cat("Using fit to convert phi into pval...\n")
    pval <- fit(prop@pairs$prop)

    cat("Correcting for multiple testing...\n")
    fdr <- p.adjust(pval, method = "BH")

    cat("Building results...\n")
    prop@pairs <- data.frame(prop@pairs, "pval" = pval, "fdr" = fdr, stringsAsFactors = FALSE)
    prop@pairs <- prop@pairs[order(prop@pairs$pval), ]
    rownames(prop@pairs) <- 1:nrow(prop@pairs)

    return(prop)

  }else{

    return(fit)
  }
}
