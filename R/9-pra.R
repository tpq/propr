#' Principal Ratio Analysis
#'
#' This function finds which feature ratios explain the most variance.
#'  This is a computationally expensive procedure that we approximate
#'  with the heuristic described below.
#'
#' This function resembles the method described by Michael Greenacre
#'  in "Variable Selection in Compositional Data Analysis Using
#'  Pairwise Logratios", except that we have modified the method
#'  to use a heuristic that scales to high-dimensional data.
#'
#' For each ratio, the heuristic will search CLR-based clusters
#'  for the best denominator, and then will search ALR-based clusters
#'  for the best numerator. It does this by dividing the
#'  transformed data into \code{nclust} clusters, calculating
#'  \code{vegan::rda} on the geometric mean of each cluster, then
#'  searching the best clusters exhaustively. The \code{ndenom}
#'  argument toggles how many best denominators to use during the
#'  next step. This process is repeated \code{ndim} times, finding
#'  that number of ratios that explain the most variance.
#'
#' @inheritParams all
#' @param ndim An integer. The number of ratios to find.
#' @param nclust An integer. The number of clusters to build from the data.
#' @param nsearch An integer. The number of clusters to search exhaustively.
#' @param ndenom An integer. The number of best denominators to use
#'  when searching for the best numerators.
#'
#' @return A list of: (1) "best", the best ratios and the variance they explain,
#'  (2) "all", all ratios tested and the variance they explain,
#'  (3) "Z", the standardized data used by the constrained PCA, and
#'  (4) "Y", the final ratios used to constrain the PCA.
#'
#' @export
pra <- function(counts, ndim = 3, nclust = 2*round(sqrt(ncol(counts))), nsearch = 3, ndenom = 4){

  packageCheck("vegan")

  searchTree <- function(data, Z, nclust = ncol(data) / 10, nsearch = 1, lrm = NULL){

    if(nsearch > nclust) stop("You cannot have more 'nsearch' than 'nclust'.")

    all0s <- apply(data, 2, function(x) all(x == 0))
    d <- dist(t(data[,!all0s]))
    h <- fastcluster::hclust(d)
    cuts <- cutree(h, nclust)

    # Spot check clades
    l1 <- sapply(1:nclust, function(cut){

      index <- names(cuts)[cuts == cut]
      clrOfCut <- rowMeans(data[, index, drop = FALSE])
      lr.try <- cbind(lrm, clrOfCut)
      rs <- vegan::rda(Z, lr.try)
      sum(rs$CCA$eig)/(sum(rs$CA$eig)+sum(rs$CCA$eig))*100
    })

    cuts.best <- order(l1, decreasing = TRUE)[1:nsearch]
    trythese <- names(cuts)[cuts %in% cuts.best]
    l2 <- sapply(trythese, function(id){

      memberOfCut <- data[,id]
      lr.try <- cbind(lrm, memberOfCut)
      rs <- vegan::rda(Z, lr.try)
      sum(rs$CCA$eig)/(sum(rs$CA$eig)+sum(rs$CCA$eig))*100
    })

    l2
  }

  # This method requires column names
  if(is.null(colnames(counts))){

    colnames(counts) <- paste0(1:ncol(counts))
  }

  # There is a maximum number of ratios that can be found
  if(ndim > min(dim(counts))-1){

    message("Alert: You have requested too many dimensions.")
    message("Alert: Retrieving all dimensions instead.")
    ndim <- min(dim(counts))-1
  }

  # Zeros will cause problems with this method
  if(any(as.matrix(counts) == 0)){

    message("Alert: Replacing 0s with next smallest value.")
    zeros <- counts == 0
    counts[zeros] <- min(counts[!zeros])
  }

  # Calculate Z used to fit vegan model
  P <- counts/sum(counts)
  rm <- apply(P,1,sum)
  cm <- apply(P,2,sum)
  Y <- as.matrix(log(P))
  mc <- t(Y) %*% as.vector(rm)
  Y <- Y - rep(1, nrow(P)) %*% t(mc)
  mr <- Y %*% as.vector(cm)
  Y <- Y - mr %*% t(rep(1, ncol(P)))
  Z <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))

  # Find k ratios that explain most variance
  bmat <- vector("list", ndim)
  lrm <- NULL # do not delete -tpq
  for(k in 1:ndim){

    # Phase I: Guess best DENOMINATOR based on CLR-transformed input
    clr <- log(counts) - rowMeans(log(counts))
    cxv <- searchTree(clr, Z, nclust, nsearch, lrm)
    topOrAll <- min(ndenom, length(cxv))
    best <- names(cxv)[order(cxv, decreasing = TRUE)][1:topOrAll]

    # Phase II: Guess best NUMERATOR based on ALR-transformed input
    bmat[[k]] <- lapply(best, function(b){

      # ALR-transform using the b-th best feature
      alr <- log(counts) - log(counts[,b])
      cxv <- searchTree(alr, Z, nclust, nsearch, lrm)
      explainedVar <- cxv[order(cxv, decreasing = TRUE)]

      # Return
      data.frame(
        "k" = k,
        "Partner" = b,
        "Pair" = names(explainedVar),
        "var" = explainedVar,
        stringsAsFactors = FALSE
      )
    })

    # Record variance explained by each ratio
    bmat[[k]] <- do.call("rbind", bmat[[k]])
    ind <- order(bmat[[k]]$var, decreasing=TRUE)
    bmat[[k]] <- bmat[[k]][ind, ]

    # Update LRM based on best ratio
    numer <- bmat[[k]][1, "Pair"]
    denom <- bmat[[k]][1, "Partner"]
    lrm <- cbind(lrm, log(counts[, numer] / counts[, denom]))

    # Update progress bar
    numTicks <- progress(k, ndim, numTicks)
  }

  # Get best ratio
  resm <- lapply(bmat, function(x) x[1, ])

  list(
    "best" = do.call("rbind", resm),
    "all" = do.call("rbind", bmat),
    "Z" = Z,
    "Y" = lrm
  )
}
