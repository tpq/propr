#' Calculate odds ratio and FDR 
#'
#' This function calls \code{\link{graflex}} for each
#' concept (i.e., column) in the database \code{K}.
#'
#' For each concept, this function calculates the odds ratio
#' and determines the false discovery rate (FDR) by counting 
#' the number of the actual OR was greater or less than a 
#' permuted OR.
#' 
#' @param A An adjacency matrix.
#' @param K A knowledge database where each row is a graph node
#'  and each column is a concept.
#' @param p An integer. The number of permutation.
#' 
#' @export
runGraflex <- function(A, K, p=100, ncores=1) {
  if (nrow(A) != nrow(K))
    stop("'A' and 'K' must have identical rows.")
  if (nrow(A) != ncol(A))
    stop("'A' must be a square matrix.")

  if (ncores == 1){

    # for each knowledge network, calculate odds ratio and FDR
    res <- lapply(1:ncol(K), function(k) {
      Gk <- K[, k] %*% t(K[, k])      # converts the k column into an adjacency matrix (genes x genes)
      graflex(A, Gk, p=p)  # this calls the graflex function implemented in Rcpp C++
    })

  }else{
    packageCheck("parallel")
    cl <- parallel::makeCluster(ncores)
    res <- parallel::parLapply(cl, 1:ncol(K), function(k) {
      Gk <- K[, k] %*% t(K[, k])
      graflex(A, Gk, p=p)
    })
    parallel::stopCluster(cl)
  }

  # parse resulting data frame
  res <- do.call("rbind", res)
  res <- cbind(res, rep(p, ncol(K)))
  res <- cbind(res, colnames(K))
  res <- as.data.frame(res)
  colnames(res) <- c("Neither", "G.only", "A.only", "Both", "Odds", "LogOR", "FDR.under", "FDR.over", "Permutes", "Concept")

  return(res)
}
