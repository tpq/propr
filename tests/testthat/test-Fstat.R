if(requireNamespace("limma", quietly = TRUE) &
   requireNamespace("SDMTools", quietly = TRUE)){

  testFmod <- function(counts, group, a = .05){

    library(propr)
    library(limma)

    M <- counts

    gr <- group
    ce=which(group == unique(group)[1])
    co=which(group == unique(group)[2])
    cere=ce
    cort=co
    nce=length(cere)
    nco=length(cort)

    design=matrix(0,dim(M)[1],2)
    design[1:length(cere),1]=rep(1,length(cere))
    design[(length(cere)+1):dim(M)[1],2]=rep(1,length(cort))

    #geometric mean "gene":
    z=exp(apply(log(M),1,mean))
    #ratios used for the hierarchical model:
    Mz=M/z

    #scale counts up by the number mean(z) to make them similar to original counts (because counts will get pseudocounted by limma)
    scaledcounts=t(Mz*mean(z))
    u=voom(scaledcounts, design=design, plot=TRUE)
    # this is only used for the parameters of the hierarchical model (the resulting count weights are incompatible with the weights taken directly from M that are used in propr)
    param=lmFit(u,design)
    param=eBayes(param)
    dz=param$df.prior
    s2z=param$s2.prior

    #now we get the usual counts and weights:
    counts=t(M[c(cere,cort),])
    v=voom(counts, design=design, plot=TRUE)
    colnames(v$weights)=rownames(M)

    #get weight of the geometric mean "gene"
    su=apply(v$weights,2,sum)
    w=t(v$weights)/su
    #weighted geometric mean
    zw=apply(w*log(M)*dim(M)[2],1,mean)
    #define weight of geometric mean:
    wz=zw/apply(log(M),1,mean)

    ###########################################################

    #Let's now assume unweighted results are in object "res", weighted results in object "Res":

    message("(0) initializing propr (unweighted):")
    res=propd(counts=M, group=gr, p = 1)

    message("(0) initializing propr (weighted):")
    Res=propd(counts=M, group=gr, p = 1, weighted = TRUE)

    #from these results, we need only theta and LRV:
    st=res@theta[,c("lrv","theta")]
    stw=Res@theta[,c("lrv","theta")]

    message("(1) unweighted moderated statistic:")

    #within-group variances wrt reference z:
    s2xz=rep(0,dim(Mz)[2])
    for (i in 1:dim(Mz)[2]){

      s1=var(log(Mz[ce,i]))
      s2=var(log(Mz[co,i]))
      s2xz[i]=((nce-1)*s1+(nco-1)*s2)/(nce+nco-1)
    }

    #moderation term for each ratio:
    mod=c(1:dim(st)[1])
    i=0
    for (j in 2:dim(Mz)[2]){
      print(j)
      for (k in 1:(j-1)){
        i=i+1
        mod[i]=2*s2z-s2xz[j]-s2xz[k]
      }
    }
    mod=mod/(st[,"lrv"]*st[,"theta"]*(1+(nce+nco)/dz))

    Fpmod=(1-st[,"theta"])/(st[,"theta"]*(1+mod))
    thetamod=1/(1+Fpmod)
    Fmod=(nce+nco-2)*Fpmod

    message("(2) weighted moderated statistic:")

    #weighted within-group variances wrt reference z:
    sw2xz=rep(0,dim(Mz)[2])
    for (i in 1:dim(Mz)[2]){

      W=v$weights[i,]*wz


      n=sum(W)
      s=sum(W^2)
      p=n-s/n

      n1=sum(W[ce])
      s1=sum(W[ce]^2)
      p1=n1-s1/n1

      n2=sum(W[co])
      s2=sum(W[co]^2)
      p2=n2-s2/n2



      sw1=SDMTools::wt.var(log(Mz[ce,i]),W[ce])
      sw2=SDMTools::wt.var(log(Mz[co,i]),W[co])
      sw2xz[i]=(p1*sw1+p2*sw2)/p

    }

    #weighted moderation term for each ratio:
    modw=c(1:dim(stw)[1])
    i=0
    for (j in 2:dim(Mz)[2]){
      print(j)
      for (k in 1:(j-1)){
        i=i+1
        modw[i]=2*s2z-sw2xz[j]-sw2xz[k]

      }
    }
    modw=modw/(stw[,"lrv"]*stw[,"theta"]*(1+(nce+nco)/dz))

    Fpmodw=(1-stw[,"theta"])/(stw[,"theta"]*(1+modw))
    thetamodw=1/(1+Fpmodw)
    Fmodw=(nce+nco-2)*Fpmodw

    message("(3) power-transformed moderated statistic:")

    #pre-calculate all x^a
    Ma=M^a

    #logratio variance function for pre-transformed data:
    lrvMa=function(x,y,a){
      N=length(x)
      s=sum((x/mean(x)-y/mean(y))^2)/((N-1)*a^2)
      return(s)
    }


    #power-transformed within-group variances wrt reference z:
    sa2xz=rep(0,dim(Ma)[2])
    for (i in 1:dim(Ma)[2]){

      sa1=lrvMa(Ma[ce,i],(z[ce])^a,a)
      sa2=lrvMa(Ma[co,i],(z[co])^a,a)

      sa2xz[i]=((nce-1)*sa1+(nco-1)*sa2)/(nce+nco-1)

    }

    #power-transformed moderation term, variance and theta for each ratio:
    moda=c(1:dim(st)[1])
    sta=matrix(0,dim(st)[1],2)
    colnames(sta)=c("lrv","theta")
    i=0
    for (j in 2:dim(Ma)[2]){
      print(j)
      for (k in 1:(j-1)){
        i=i+1

        moda[i]=2*s2z-sa2xz[j]-sa2xz[k]
        sxy1=lrvMa(Ma[ce,j],Ma[ce,k],a)
        sxy2=lrvMa(Ma[co,j],Ma[co,k],a)
        sta[i,"lrv"]=lrvMa(Ma[,j],Ma[,k],a)
        sta[i,"theta"]=((nce-1)*sxy1+(nco-1)*sxy2)/((nce+nco-1)*sta[i,"lrv"])

      }
    }
    moda=moda/(sta[,"lrv"]*sta[,"theta"]*(1+(nce+nco)/dz))
    Fpmoda=(1-sta[,"theta"])/(sta[,"theta"]*(1+moda))
    thetamoda=1/(1+Fpmoda)
    Fmoda=(nce+nco-2)*Fpmoda

    message("(4) weighted power-transformed moderated statistic:")

    #pre-calculate all x^a
    Ma=M^a

    #weighted logratio variance function for pre-transformed data:
    #note that this returns a precision-weighted variance like in the SDMTools package:
    lrvMaw=function(x,y,a,W){

      n=sum(W)
      s=sum(W^2)
      p=n-s/n

      N=length(x)
      w=W/n

      s=sum(W*(x/(N*mean(w*x))-y/(N*mean(w*y)))^2)/(p*a^2)
      return(s)
    }

    #weighted power-transformed within-group variances wrt reference z:
    swa2xz=rep(0,dim(Ma)[2])
    for (i in 1:dim(Ma)[2]){

      W=v$weights[i,]*wz
      swa1=lrvMaw(Ma[ce,i],(z[ce])^a,a,W[ce])
      swa2=lrvMaw(Ma[co,i],(z[co])^a,a,W[co])

      n=sum(W)
      s=sum(W^2)
      p=n-s/n

      n1=sum(W[ce])
      s1=sum(W[ce]^2)
      p1=n1-s1/n1

      n2=sum(W[co])
      s2=sum(W[co]^2)
      p2=n2-s2/n2

      swa2xz[i]=(p1*swa1+p2*swa2)/p

    }

    #weighted and power-transformed moderation term, variance and theta for each ratio:
    modwa=c(1:dim(st)[1])
    stwa=matrix(0,dim(st)[1],2)
    colnames(stwa)=c("lrv","theta")
    i=0
    for (j in 2:dim(Ma)[2]){
      print(j)
      for (k in 1:(j-1)){
        i=i+1

        W=v$weights[j,]*v$weights[k,]

        n=sum(W)
        s=sum(W^2)
        p=n-s/n

        n1=sum(W[ce])
        s1=sum(W[ce]^2)
        p1=n1-s1/n1

        n2=sum(W[co])
        s2=sum(W[co]^2)
        p2=n2-s2/n2

        modwa[i]=2*s2z-swa2xz[j]-swa2xz[k]
        swxy1=lrvMaw(Ma[ce,j],Ma[ce,k],a,W[ce])
        swxy2=lrvMaw(Ma[co,j],Ma[co,k],a,W[co])
        stwa[i,"lrv"]=lrvMaw(Ma[,j],Ma[,k],a,W)
        stwa[i,"theta"]=(p1*swxy1+p2*swxy2)/(p*stwa[i,"lrv"])

      }
    }
    modwa=modwa/(stwa[,"lrv"]*stwa[,"theta"]*(1+(nce+nco)/dz))
    Fpmodwa=(1-stwa[,"theta"])/(stwa[,"theta"]*(1+modwa))
    thetamodwa=1/(1+Fpmodwa)
    Fmodwa=(nce+nco-2)*Fpmodwa

    return(list(thetamod, thetamodw, thetamoda, thetamodwa,
                Fmod, Fmodw, Fmoda, Fmodwa))
  }

  library(propr)

  data(iris)
  keep <- iris$Species %in% c("setosa", "versicolor")
  counts <- iris[keep, 1:4] * 10
  group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

  pd.nn <- propd(counts, group)
  pd.wn <- propd(counts, group, weighted = TRUE)
  pd.na <- propd(counts, group, alpha = .05)
  pd.wa <- propd(counts, group, weighted = TRUE, alpha = .05)

  pd.nn <- updateF(pd.nn, moderated = TRUE)
  pd.wn <- updateF(pd.wn, moderated = TRUE)
  pd.na <- updateF(pd.na, moderated = TRUE)
  pd.wa <- updateF(pd.wa, moderated = TRUE)

  ref <- testFmod(counts, group, a = .05)

  test_that("updateF matches code provided by Ionas", {

    expect_equal(
      pd.nn@theta$theta_mod,
      ref[[1]]
    )

    expect_equal(
      pd.nn@theta$Fstat,
      ref[[5]]
    )

    expect_equal(
      pd.wn@theta$theta_mod,
      ref[[2]]
    )

    expect_equal(
      pd.wn@theta$Fstat,
      ref[[6]]
    )

    expect_equal(
      pd.na@theta$theta_mod,
      ref[[3]]
    )

    expect_equal(
      pd.na@theta$Fstat,
      ref[[7]]
    )

    expect_equal(
      pd.wa@theta$theta_mod,
      ref[[4]]
    )

    expect_equal(
      pd.wa@theta$Fstat,
      ref[[8]]
    )
  })
}
