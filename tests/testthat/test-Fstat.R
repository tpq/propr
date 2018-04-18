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

    scaledcounts=t(Mz*mean(z))
    u=voom(scaledcounts, design=design, plot=TRUE)
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
    zw=exp(apply(w*log(M)*dim(M)[2],1,mean))
    #define weight of geometric mean:
    #wz=zw/z
    wz=log(zw)/log(z)

    res=propd(counts=M, group=gr, p = 1)
    Res=propd(counts=M, group=gr, p = 1, weighted = TRUE)

    #from these results, we need only theta and LRV:
    st=res@theta[,c("lrv","theta")]
    stw=Res@theta[,c("lrv","theta")]

    #print("(1) unweighted moderated statistic")

    mod=dz*s2z/st[,"lrv"]
    Fpmod=(1-st[,"theta"])*(dz+nce+nco)/((nce+nco)*st[,"theta"]+mod)
    thetamod=1/(1+Fpmod)
    Fmod=(nce+nco+dz-2)*Fpmod

    #print("(2) weighted moderated statistic")

    modw=dz*s2z/stw[,"lrv"]
    Fpmodw=(1-stw[,"theta"])*(dz+nce+nco)/((nce+nco)*stw[,"theta"]+modw)
    thetamodw=1/(1+Fpmodw)
    Fmodw=(nce+nco+dz-2)*Fpmodw

    #print("(3) power-transformed moderated statistic")

    #logratio variance function for pre-transformed data:
    Ma=M^a
    lrvMa=function(x,y,a){
      N=length(x)
      s=sum((x/mean(x)-y/mean(y))^2)/((N-1)*a^2)
      return(s)
    }

    sta=matrix(0,dim(st)[1],2)
    colnames(sta)=c("lrv","theta")
    i=0
    for (j in 2:dim(Ma)[2]){
      for (k in 1:(j-1)){
        i=i+1
        sxy1=lrvMa(Ma[ce,j],Ma[ce,k],a)
        sxy2=lrvMa(Ma[co,j],Ma[co,k],a)
        sta[i,"lrv"]=lrvMa(Ma[,j],Ma[,k],a)
        sta[i,"theta"]=((nce-1)*sxy1+(nco-1)*sxy2)/((nce+nco-1)*sta[i,"lrv"])
      }
    }

    moda=dz*s2z/sta[,"lrv"]
    Fpmoda=(1-sta[,"theta"])*(dz+nce+nco)/((nce+nco)*sta[,"theta"]+moda)
    thetamoda=1/(1+Fpmoda)
    Fmoda=(nce+nco+dz-2)*Fpmoda

    #print("(4) weighted power-transformed moderated statistic")

    #weighted logratio variance function for pre-transformed data:
    #note that this returns a precision-weighted variance like in the SDMTools package:
    Ma=M^a
    lrvMaw=function(x,y,a,W){
      n=sum(W)
      s=sum(W^2)
      p=n-s/n
      N=length(x)
      w=W/n
      s=sum(W*(x/(N*mean(w*x))-y/(N*mean(w*y)))^2)/(p*a^2)
      return(s)
    }

    stwa=matrix(0,dim(st)[1],2)
    colnames(stwa)=c("lrv","theta")
    i=0
    for (j in 2:dim(Ma)[2]){
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
        swxy1=lrvMaw(Ma[ce,j],Ma[ce,k],a,W[ce])
        swxy2=lrvMaw(Ma[co,j],Ma[co,k],a,W[co])
        stwa[i,"lrv"]=lrvMaw(Ma[,j],Ma[,k],a,W)
        stwa[i,"theta"]=(p1*swxy1+p2*swxy2)/(p*stwa[i,"lrv"])
      }
    }

    modwa=dz*s2z/stwa[,"lrv"]
    Fpmodwa=(1-stwa[,"theta"])*(dz+nce+nco)/((nce+nco)*stwa[,"theta"]+modwa)
    thetamodwa=1/(1+Fpmodwa)
    Fmodwa=(nce+nco+dz-2)*Fpmodwa

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
