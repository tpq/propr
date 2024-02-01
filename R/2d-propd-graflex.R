#' Permute Odds Ratio
#'
#' This function permutes \code{p} odds ratios for the
#'  edge overlap between two non-random graphs.
#'  It does this by randomly shuffling the rows and
#'  columns of \code{A} (jointly, thus preserving
#'  the degree distribution).
#'
#' Note that this function calculates overlap for the
#'  lower-left triangle of the input matrices.
#' @param A,G An adjacency matrix.
#' @param p An integer. The number of overlaps to permute.
permuteOR <- function(A, G, p = 500) {
  Gstar <- G[lower.tri(G)]
  res <- lapply(1:p, function(i) {
    # Shuffle the adjacency matrix
    index <- sample(1:ncol(A))
    A <- A[index, index]
    Astar <- A[lower.tri(A)]
    getOR(Astar, Gstar)
  })

  do.call("rbind", res)
}

#' Tabulate Overlap
#'
#' This function tabulates the overlap between
#'  two vectors or two adjacency matrices.
#'  It is a faster version of \code{table} that
#'  only supports binary input.
#' @param A,G A vector or adjacency matrix.
#' @return A table of overlap.
binTab <- function(A, G) {
  diff <- A != G
  only1 <- A[diff]
  b <- sum(only1)
  c <- length(only1) - b

  same <- !diff
  double1 <- A[same]
  a <- sum(double1)
  d <- length(double1) - a

  matrix(c(d, b, c, a), 2, 2)
}

#' Calculate Odds Ratio
#'
#' This function calculates the overlap between
#'  two vectors or two adjacency matrices.
#'  It returns the OR as well as other metrics.
#' @inheritParams binTab
#' @return A \code{data.frame} of results.
getOR <- function(A, G) {
  tab <- binTab(A, G)
  or <- (tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1])
  data.frame(
    "Neither" = tab[1, 1],
    "G.only" = tab[1, 2],
    "A.only" = tab[2, 1],
    "Both" = tab[2, 2],
    "Odds" = or,
    "LogOR" = log(or)
  )
}

#' Calculate Odds Ratio
#'
#' This function calculates an odds ratio for the
#'  edge overlap between two non-random graphs.
#'
#' Note that this function calculates overlap for the
#'  lower-left triangle of the input matrices.
#' @inheritParams permuteOR
calculateOR <- function(A, G) {
  Astar <- A[lower.tri(A)]
  Gstar <- G[lower.tri(G)]
  getOR(Astar, Gstar)
}

#' Calculate Odds Ratio FDR
#'
#' This function calculates the false discovery rate (FDR)
#'  for over- and under-enrichment by counting the number of
#'  times the actual OR was greater than
#'  (or less than) a permuted OR.
#' @param actual A result from \code{\link{calculateOR}}.
#' @param permuted A result from \code{\link{permuteOR}}.
#' @return A \code{data.frame} of the FDRs for over-
#'  and under- enrichment.
getFDR <- function(actual, permuted) {
  actual$FDR.under <-
    sum(permuted$Odds <= actual$Odds) / nrow(permuted)
  actual$FDR.over <-
    sum(permuted$Odds >= actual$Odds) / nrow(permuted)
  actual
}

#' Permute FDR for Multiple Concepts
#'
#' This function calls \code{\link{permuteOR}} for each
#'  concept (i.e., column) in the database \code{K}.
#'
#' For each concept, this function calculates the
#'  false discovery rate (FDR) by counting the number of
#'  times the actual OR was greater than
#'  (or less than) a permuted OR.
#' @inheritParams permuteOR
#' @param A An adjacency matrix.
#' @param K A knowledge database where each row is a graph node
#'  and each column is a concept.
#' @export
runGraflex <- function(A, K, p = 500) {
  if (nrow(A) != nrow(K))
    stop("'A' and 'K' must have identical rows.")

  numTicks <- 0
  res <- lapply(1:ncol(K), function(k) {
    numTicks <<- progress(k, ncol(K), numTicks)
    Gk <- K[, k] %*% t(K[, k])
    actual <- calculateOR(A, Gk)
    permuted <- permuteOR(A, Gk, p = p)
    actual <- getFDR(actual, permuted)
    actual$Permutes <- p
    actual$Concept <- colnames(K)[k]
    actual
  })

  do.call("rbind", res)
}
