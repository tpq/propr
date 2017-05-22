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

  lrvMa=function(x,y,a){
    N=length(x)
    s=sum((x/mean(x)-y/mean(y))^2)/(a^2*(N-1))
    return(s)
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
      ast[i,4]=lrvMa(Ma[,j],Ma[,k],a)
      ast[i,5]=lrvMa(Ma[cere,j],Ma[cere,k],a)
      ast[i,6]=lrvMa(Ma[cort,j],Ma[cort,k],a)
      ast[i,7]=mean(Ma[cere,j]-Ma[cere,k])/a
      ast[i,8]=mean(Ma[cort,j]-Ma[cort,k])/a
      ast[i,3]=((n1-1)*ast[i,5]+(n2-1)*ast[i,6])/((n1+n2-1)*ast[i,4])
      i=i+1
    }
  }

  return(as.data.frame(ast))
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

#' Calculate alpha Theta
#'
#' Calculate differential proportionality measure, theta,
#'  using the Box-Cox transformation method.
#'
#' @inheritParams propd
alphaTheta <- function(counts, group, alpha){

  ct <- as.matrix(counts)
  lrv <- boxRcpp(ct[], alpha)

  if(length(unique(group)) != 2) stop("Please use two groups.")
  if(length(group) != nrow(counts)) stop("Too many or too few group labels.")
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]

  lrv1 <- boxRcpp(ct[group1,], alpha)
  lrv2 <- boxRcpp(ct[group2,], alpha)
  n1 <- sum(group1)
  n2 <- sum(group2)

  theta <- ((n1-1) * lrv1 + (n2-1) * lrv2) / ((n1+n2-1) * lrv)

  labels <- labRcpp(ncol(counts))
  return(
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      "theta" = theta,
      "lrv" = lrv,
      "lrv1" = lrv1,
      "lrv2" = lrv2
    ))
}
