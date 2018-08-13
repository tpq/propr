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
#' @param X A data.frame or matrix. A "count matrix" with
#'  subjects as rows and features as columns.
#' @param check A logical. If TRUE, function first checks
#'  for negative and NA values.
#' @return Returns a matrix containing the variance of the
#'  log of the ratios.
proprVLR <- function(X, check = FALSE){

  if(check){

    if(any(X <= 0)) stop("non-positive values found")
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
#' @return A matrix. Returns the centered log-ratio
#'  transformation of \code{X}.
proprCLR <- function(X, check = FALSE){

  if(check){

    if(any(X <= 0)) stop("non-positive values found")
    if(any(is.na(X))) stop("NA values found")
  }

  logX <- log(X)
  return(sweep(logX, 1, rowMeans(logX), "-")) # subtract out the means
}

#' Calculates the additive log-ratio transformation.
#'
#' Provided for backend use.
#'
#' @param ivar A numeric scalar. Specifies feature to use
#'  as reference for additive log-ratio transformation.
#' @inheritParams proprVLR
#' @return A matrix. Returns the additive log-ratio
#'  transformation of \code{X}.
proprALR <- function(X, ivar, check = FALSE){

  if(check){

    if(any(X <= 0)) stop("non-positive values found")
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
#' @return A vector. Returns the lower left triangle of a
#'  proportionality matrix.
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

#' Calculate Theta
#'
#' Do not use this function. For testing purposes only.
#'
#' @inheritParams propd
calculateTheta_old <- function(counts, group){

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  cere=rownames(counts)[group1]
  cort=rownames(counts)[group2]
  n1=length(cere)
  n2=length(cort)

  st=matrix(0,(dim(counts)[2]*(dim(counts)[2]-1)/2),8)
  colnames(st)=c("gene1","gene2","theta","lrv","lrv1","lrv2","lrm1","lrm2")
  i=1
  for(j in 2:dim(counts)[2]){
    for(k in 1:(j-1)){

      st[i,1]=as.integer(j)
      st[i,2]=as.integer(k)
      st[i,4]=var(log(counts[,j]/counts[,k]))
      st[i,5]=var(log(counts[cere,j]/counts[cere,k]))
      st[i,6]=var(log(counts[cort,j]/counts[cort,k]))
      st[i,7]=mean(log(counts[cere,j]/counts[cere,k]))
      st[i,8]=mean(log(counts[cort,j]/counts[cort,k]))
      st[i,3]=((n1-1)*st[i,5]+(n2-1)*st[i,6])/((n1+n2-1)*st[i,4])
      i=i+1
    }
  }

  return(as.data.frame(st))
}

#' Calculate Weighted Theta
#'
#' Do not use this function. For testing purposes only.
#'
#' @inheritParams propd
calculateThetaW_old <- function(counts, group){

  packageCheck("limma")

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  cere=rownames(counts)[group1]
  cort=rownames(counts)[group2]
  n1=length(cere)
  n2=length(cort)

  design=matrix(0,dim(counts)[1],2)
  design[1:length(cere),1]=rep(1,length(cere))
  design[(length(cere)+1):dim(counts)[1],2]=rep(1,length(cort))
  v=limma::voom(t(counts), design=design, plot=TRUE)
  w=t(v$weights)

  st=matrix(0,(dim(counts)[2]*(dim(counts)[2]-1)/2),8)
  colnames(st)=c("gene1","gene2","theta","lrv","lrv1","lrv2","lrm1","lrm2")
  i=1
  for(j in 2:dim(counts)[2]){
    for(k in 1:(j-1)){

      W=w[,j]*w[,k]
      n=sum(W)
      s=sum(W^2)
      p=n-s/n

      n1=sum(W[group1])
      s1=sum(W[group1]^2)
      p1=n1-s1/n1

      n2=sum(W[group2])
      s2=sum(W[group2]^2)
      p2=n2-s2/n2

      st[i,1]=as.integer(j)
      st[i,2]=as.integer(k)
      st[i,4]=wtvRcpp(log(counts[,j]/counts[,k]), W)
      st[i,5]=wtvRcpp(log(counts[cere,j]/counts[cere,k]), W[group1])
      st[i,6]=wtvRcpp(log(counts[cort,j]/counts[cort,k]), W[group2])
      st[i,7]=wtmRcpp(log(counts[cere,j]/counts[cere,k]), W[group1])
      st[i,8]=wtmRcpp(log(counts[cort,j]/counts[cort,k]), W[group2])
      st[i,3]=(p1*st[i,5]+p2*st[i,6])/(p*st[i,4])
      i=i+1
    }
  }

  return(as.data.frame(st))
}

#' Permute Theta
#'
#' Do not use this function. For testing purposes only.
#'
#' @inheritParams propd
permuteTheta_old <- function(counts, group, p){

  theta <- calculateTheta_old(counts, group)

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  cere=rownames(counts)[group1]
  cort=rownames(counts)[group2]
  n1=length(cere)
  n2=length(cort)

  ptheta=matrix(0,dim(theta)[1],p)
  for(k in 1:p){

    message("Permuting theta for iteration: ", k)
    for(i in 1:dim(theta)[1]){

      # Sample group1 vector to mimic pairmutate
      index <- 1:(n1+n2)
      grp <- sample(group1)
      perm <- c(index[grp], index[!grp])

      # Continue with original code
      st5=var(log(counts[perm[1:n1],theta[i,1]]/counts[perm[1:n1],theta[i,2]]))
      st6=var(log(counts[perm[(n1+1):(n1+n2)],theta[i,1]]/counts[perm[(n1+1):(n1+n2)],theta[i,2]]))
      ptheta[i,k]=((n1-1)*st5+(n2-1)*st6)/((n1+n2-1)*theta[i,4])
    }
  }

  return(ptheta)
}

#' Calculate alpha Theta
#'
#' Do not use this function. For testing purposes only.
#'
#' @inheritParams propd
alphaTheta_old <- function(counts, group, alpha){

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  cere=rownames(counts)[group1]
  cort=rownames(counts)[group2]
  n1=length(cere)
  n2=length(cort)

  lrva=function(group,a,x,y){

    N=length(group)
    lrv=sum(((x[group]-mean(x[group]))/mean(x)-(y[group]-mean(y[group]))/mean(y))^2)
    lrv=lrv/((N-1)*a^2)
    return(lrv)
  }

  lrma=function(group,a,x,y){

    m=mean(x[group]/mean(x)- y[group]/mean(y))
    gr2=setdiff(c(1:length(x)),group)
    C=mean(x[group]-y[group])+mean(x[gr2]-y[gr2])
    lrm=(C/2+m)/a
    return(lrm)
  }

  a=alpha
  Ma=as.data.frame(counts^a) # pre-calculate all x^a
  ast=matrix(0,(dim(Ma)[2]*(dim(Ma)[2]-1)/2),8)
  colnames(ast)=c("gene1","gene2","atheta","alrv","alrv1","alrv2","alrm1","alrm2")
  i=1
  for (j in 2:dim(Ma)[2]){
    for (k in 1:(j-1)){

      ast[i,1]=as.integer(j)
      ast[i,2]=as.integer(k)
      ast[i,4]=lrva(which(group1 | group2), a, Ma[,j], Ma[,k])
      ast[i,5]=lrva(which(group1), a, Ma[,j], Ma[,k])
      ast[i,6]=lrva(which(group2), a, Ma[,j], Ma[,k])
      ast[i,7]=lrma(which(group1), a, Ma[,j], Ma[,k])
      ast[i,8]=lrma(which(group2), a, Ma[,j], Ma[,k])
      ast[i,3]=((n1-1)*ast[i,5]+(n2-1)*ast[i,6])/((n1+n2-1)*ast[i,4])
      i=i+1
    }
  }

  return(as.data.frame(ast))
}

#' Calculate Weighted Theta
#'
#' Do not use this function. For testing purposes only.
#'
#' @param weights A matrix. Pre-computed \code{limma}-based
#'  weights. Optional parameter.
#' @inheritParams propd
alphaThetaW_old <- function(counts, group, alpha, weights){

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  cere=as.numeric(rownames(counts)[group1])
  cort=as.numeric(rownames(counts)[group2])
  n1=length(cere)
  n2=length(cort)

  lrvaw=function(group,a,x,y,W){

    n=sum(W[group])
    s=sum(W[group]^2)
    p=n-s/n
    w=W[group]/n

    lrv=sum(W[group]*((x[group]-sum(w*x[group]))/(sum(W*x)/sum(W))-(y[group]-sum(w*y[group]))/(sum(W*y)/sum(W)))^2)
    lrv=lrv/(p*a^2)
    return(lrv)
  }

  lrmaw=function(group,a,x,y,W){

    w=W[group]/sum(W[group])
    m=sum(w*(x[group]/(sum(W*x)/sum(W))-y[group]/(sum(W*y)/sum(W))))

    gr2=setdiff(c(1:length(x)),group)
    w2=W[gr2]/sum(W[gr2])
    C=sum(w*(x[group]-y[group]))+sum(w2*(x[gr2]-y[gr2]))

    lrm=(C/2+m)/a
    return(lrm)
  }

  a=alpha
  Ma=counts^a
  llt=ncol(counts)*(ncol(counts)-1)/2
  modwa=c(1:llt)
  stwa=matrix(0,llt,6)
  colnames(stwa)=c("lrv","theta", "awlrv1", "awlrv2", "awlrm1", "awlrm2")
  i=0
  for (j in 2:dim(Ma)[2]){
    for (k in 1:(j-1)){

      i=i+1
      W=weights[,j]*weights[,k]
      n=sum(W)
      s=sum(W^2)
      p=n-s/n
      n1=sum(W[cere])
      s1=sum(W[cere]^2)
      p1=n1-s1/n1
      n2=sum(W[cort])
      s2=sum(W[cort]^2)
      p2=n2-s2/n2
      swxy1=lrvaw(which(group1), a, Ma[,j], Ma[,k], W)
      swxy2=lrvaw(which(group2), a, Ma[,j], Ma[,k], W)
      stwa[i,"lrv"]=lrvaw(which(group1 | group2), a, Ma[,j], Ma[,k], W)
      stwa[i,"theta"]=(p1*swxy1+p2*swxy2)/(p*stwa[i,"lrv"])
      stwa[i,"awlrv1"] = swxy1
      stwa[i,"awlrv2"] = swxy2
      stwa[i,"awlrm1"]=lrmaw(which(group1), a, Ma[,j], Ma[,k], W)
      stwa[i,"awlrm2"]=lrmaw(which(group2), a, Ma[,j], Ma[,k], W)
    }
  }

  return(as.data.frame(stwa))
}

#' Calculate FDR Cutoffs
#'
#' Do not use this function. For testing purposes only.
#'
#' @param theta A \code{data.frame}. The result of a
#'  \code{calculateTheta} call.
#' @param ptheta A \code{data.frame}. The result of a
#'  \code{permuteTheta} call.
#' @inheritParams propd
calculateFDR_old <- function(theta, ptheta, cutoff = seq(.6, .9, .1)){

  FDR=matrix(0,length(cutoff),2)
  colnames(FDR)=c("cutoff","FDR")
  FDR[,1]=cutoff
  for(i in 1:dim(FDR)[1]){
    true=length(which(theta$theta<FDR[i,"cutoff"]))
    rand=sum(apply(ptheta,1,function(x){length(which(x<FDR[i,"cutoff"]))}))/ncol(ptheta)
    FDR[i,2]=rand/true
  }

  return(FDR)
}

#' Permute Theta
#'
#' Permute differential proportionality measure, theta.
#'
#' This function randomizes group membership \code{p*nfeat} times.
#'
#' @inheritParams propd
permuteTheta <- function(counts, group, p = 64){

  ct <- as.matrix(counts)
  lrv <- lltRcpp(vlrRcpp(ct[]))

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few groups.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]
  n1 <- sum(group1)
  n2 <- sum(group2)

  ptheta <- as.data.frame(matrix(0, ncol(counts)*(ncol(counts)-1)/2, p))

  for(i in 1:p){

    # Generate pairwise permutations in C++
    # Returns half-matrix vectors containing lrv1 and lrv2
    message("Permuting theta for iteration: ", i)
    lrv12 <- pairmutate(ct, group1)
    lrv1 <- lrv12[[1]]
    lrv2 <- lrv12[[2]]

    # Use lrv1 and lrv2 vectors to calculate theta
    ptheta[, i] <- ((n1-1) * lrv1 + (n2-1) * lrv2) / ((n1+n2-1) * lrv)
  }

  return(ptheta)
}

#' Permute Theta
#'
#' Permute differential proportionality measure, theta.
#'
#' This function randomizes all feature vectors \code{p} times.
#'
#' @inheritParams propd
permuteTheta_naive <- function(counts, group, p = 64){

  ct <- as.matrix(counts)
  ptheta <- as.data.frame(matrix(0, ncol(counts)*(ncol(counts)-1)/2, p))

  for(i in 1:p){

    message("Permuting theta for iteration: ", i)
    ct <- apply(ct, 2, sample)
    ptheta[, i] <- calculateTheta(ct, group)$theta
  }

  return(ptheta)
}

#' Permute Theta
#'
#' Permute differential proportionality measure, theta.
#'
#' This function randomizes group labels \code{p} times.
#'
#' @inheritParams propd
permuteTheta_prime <- function(counts, group, p = 64){

  ct <- as.matrix(counts)
  ptheta <- as.data.frame(matrix(0, ncol(counts)*(ncol(counts)-1)/2, p))

  for(i in 1:p){

    message("Permuting theta for iteration: ", i)
    ct.i <- ct[sample(1:nrow(ct)), ]
    ptheta[, i] <- calculateTheta(ct.i, group)$theta
  }

  return(ptheta)
}

#' Permute Theta
#'
#' Permute differential proportionality measure, theta.
#'
#' For back-end use only.
#'
#' @inheritParams propd
permuteTheta_false <- function(counts, group, p = 64){

  ct <- as.matrix(counts)
  ptheta <- as.data.frame(matrix(0, ncol(counts)*(ncol(counts)-1)/2, p))

  for(i in 1:p){

    message("Permuting theta for iteration: ", i)
    ptheta[, i] <- calculateTheta(ct, group)$theta
  }

  return(ptheta)
}

#' Calculate FDR Cutoffs
#'
#' This function uses the result of a \code{permuteTheta} call
#'  to calculate false discovery rate (FDR) cutoffs.
#'
#' @param theta A \code{data.frame}. The result of a
#'  \code{calculateTheta} call.
#' @param ptheta A \code{data.frame}. The result of a
#'  \code{permuteTheta} call.
#' @inheritParams propd
calculateFDR <- function(theta, ptheta, cutoff = seq(.6, .9, .1)){

  out <- as.data.frame(matrix(0, length(cutoff), 2))
  colnames(out) <- c("cutoff", "FDR")
  out$cutoff <- cutoff

  for(i in 1:nrow(out)){

    message("Calculating FDR for cutoff: ", out$cutoff[i])
    true <- sum(theta$theta < out[i, "cutoff"])
    rand <- mean(apply(ptheta, 2, function(x) sum(x < out[i, "cutoff"])))
    out[i, "FDR"] <- rand / true
  }

  return(out)
}
