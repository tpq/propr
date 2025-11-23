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
runGraflex <- function(A, K, p = 100, ncores = 1) {
  NVTX_PUSH("runGraflex", 0)

  # basic checks
  if (nrow(A) != nrow(K)) {
    NVTX_POP()  # runGraflex
    stop("'A' and 'K' must have identical rows.")
  }
  if (nrow(A) != ncol(A)) {
    NVTX_POP()  # runGraflex
    stop("'A' must be a square matrix.")
  }
  if (all(rownames(A) != rownames(K))) {
    NVTX_POP()  # runGraflex
    stop("'A' and 'K' must have the same row names.")
  }

  NVTX_PUSH("compute_graflex", 0)
  if (ncores == 1) {
    NVTX_PUSH("serial_graflex", 0)
    # for each knowledge network, calculate odds ratio and FDR
    res <- lapply(1:ncol(K), function(k) {
      graflex(A, K[, k], p = p)  # calls the modified graflex function implemented in Rcpp C++
    })
    NVTX_POP()  # serial_graflex
  } else {
    NVTX_PUSH("parallel_graflex", 0)
    packageCheck("parallel")

    NVTX_PUSH("makeCluster", 0)
    cl <- parallel::makeCluster(ncores)
    NVTX_POP()  # makeCluster

    NVTX_PUSH("clusterExport", 0)
    parallel::clusterExport(cl, varlist = c("A", "K", "p"), envir = environment())
    NVTX_POP()  # clusterExport

    NVTX_PUSH("parLapply", 0)
    res <- parallel::parLapply(cl, 1:ncol(K), function(k) {
      graflex(A, K[, k], p = p)
    })
    NVTX_POP()  # parLapply

    NVTX_PUSH("stopCluster", 0)
    parallel::stopCluster(cl)
    NVTX_POP()  # stopCluster

    NVTX_POP()  # parallel_graflex
  }
  NVTX_POP()  # compute_graflex

  NVTX_PUSH("parse_results", 0)
  # parse resulting data frame
  res <- do.call("rbind", res)
  res <- cbind(res, rep(p, ncol(K)))
  res <- cbind(res, colnames(K))
  res <- as.data.frame(res)
  colnames(res) <- c(
    "Neither", "G.only", "A.only", "Both",
    "Odds", "LogOR", "FDR.under", "FDR.over",
    "Permutes", "Concept"
  )
  # change the values to numeric, except for the concept column
  res[, 1:9] <- lapply(res[, 1:9], as.numeric)
  NVTX_POP()  # parse_results

  NVTX_POP()  # runGraflex
  return(res)
}
