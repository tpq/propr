#' Calculate proportionality metric phi (Lovell 2015).
#'
#' Provided for backend use.
#'
#' @inheritParams phit
#' @return Returns proportionality matrix.
proprPhit <- function(counts, symmetrize = TRUE){

  # Replace zeroes with next smallest number
  counts[counts == 0] <- unique(sort(as.matrix(counts)))[2]

  # Calculate the variance of the log-ratio (i.e., "variation array")
  mat <- proprVLR(counts)
  colnames(mat) <- colnames(counts)
  rownames(mat) <- colnames(counts)

  # Perform clr-transformation of "count matrix"
  counts.clr <- proprCLR(counts)

  # Sweep out feature clr variance from the variation array
  for(i in 1:ncol(mat)){

    mat[, i] <- mat[, i] / stats::var(counts.clr[, i])
  }

  # Symmetrize matrix if symmetrize = TRUE
  if(symmetrize) mat <- proprSym(mat)

  return(mat)
}

#' Calculate proportionality metric rho (Erb 2016).
#'
#' Provided for backend use.
#'
#' @inheritParams phit
#' @inheritParams perb
#' @return Returns proportionality matrix.
proprPerb <- function(counts, ivar = 0){

  # Replace zeroes with next smallest number
  counts[counts == 0] <- unique(sort(as.matrix(counts)))[2]

  # Calculate the variance of the log-ratio (i.e., "variation array")
  mat <- proprVLR(counts)
  colnames(mat) <- colnames(counts)
  rownames(mat) <- colnames(counts)

  # Perform *lr-transformation of "count matrix"
  if(ivar != 0){

    mat <- mat[-ivar, -ivar] # returns one less dimension
    counts.lr <- proprALR(counts, ivar = ivar) # returns one less dimension

  }else{

    counts.lr <- proprCLR(counts)
  }

  # Sweep out feature *lr variance from the variation array
  var.lr <- apply(counts.lr, 2, stats::var)
  for(i in 1:ncol(mat)){

    for(j in 1:nrow(mat)){

      # Calculate: rho = 1 - (var(x - y))/(var(x) + var(y))
      mat[i, j] <- 1 - (mat[i, j] / (var.lr[i] + var.lr[j]))
    }
  }

  return(mat)
}

#' Calculates the variance of the log of the ratios.
#'
#' Provided for backend use.
#'
#' @param X A data.frame or matrix. A "count matrix" with subjects as rows and features as columns.
#' @param check A logical. If TRUE, function first checks for negative and NA values.
#' @return Returns a matrix containing the variance of the log of the ratios.
proprVLR <- function(X, check = FALSE){

  if(check){

    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }

  X <- log(X)
  X <- stats::var(X)
  X.diag <- diag(X)

  for(col in 1:ncol(X)){

    X[, col] <- -2 * X[, col] + X.diag + X.diag[col]
  }

  return(X)
}

#' Calculates the centered log-ratio transformation.
#'
#' Provided for backend use.
#'
#' @inheritParams proprVLR
#' @return A matrix. Returns the centered log-ratio transformation of \code{X}.
proprCLR <- function(X, check = FALSE){

  if(check){

    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }

  logX <- log(X)
  return(sweep(logX, 1, rowMeans(logX), "-")) # subtract out the means
}

#' Calculates the additive log-ratio transformation.
#'
#' Provided for backend use.
#'
#' @param ivar A numeric scalar. Specificies feature to use as reference for additive log-ratio transformation.
#' @inheritParams proprVLR
#' @return A matrix. Returns the additive log-ratio transformation of \code{X}.
proprALR <- function(X, ivar, check = FALSE){

  if(check){

    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }

  logX <- log(X[, -ivar])
  return(sweep(logX, 1, log(X[, ivar]), "-")) # subtract out the ivar
}

#' Recasts proportionality matrix as a table of feature pairs.
#'
#' Provided for backend use.
#'
#' @param prop A data.frame or matrix. A proportionality matrix.
#' @return A data.frame. Returns a table of feature pairs.
proprPairs <- function(prop){

  if(identical(dim(prop), as.integer(c(1, 1)))){

    return(data.frame())
  }

  index.i <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  index.j <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  index.prop <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  counter <- 1

  for(j in 2:nrow(prop)){

    for(i in 1:(j-1)){

      index.i[counter] <- i
      index.j[counter] <- j
      index.prop[counter] <- prop[j, i]
      counter <- counter + 1
    }
  }

  result <- data.frame("feature1" = rownames(prop)[index.i],
                       "feature2" = rownames(prop)[index.j],
                       "prop" = index.prop,
                       stringsAsFactors = FALSE)

  final <- result[order(result$prop), ]
  rownames(final) <- 1:nrow(final)

  return(final)
}

#' Retrieve the lower left triangle of a proportionality matrix.
#'
#' Provided for backend use.
#'
#' @inheritParams proprPairs
#' @return A vector. Returns the lower left triangle of a proportionality matrix.
proprTri <- function(prop){

  result <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  counter <- 1

  for(j in 2:nrow(prop)){

    for(i in 1:(j-1)){

      result[counter] <- prop[j, i]
      counter <- counter + 1
    }
  }

  return(result)
}

#' Symmetrizes a proportionality matrix.
#'
#' Provided for backend use.
#'
#' @inheritParams proprPairs
#' @return A matrix. Returns a symmetrized proportionality matrix.
proprSym <- function(prop){

  for(j in 2:nrow(prop)){

    for(i in 1:(j-1)){

      prop[i, j] <- prop[j, i]
    }
  }

  return(prop)
}
