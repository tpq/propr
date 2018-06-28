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
    resa=propd(counts=M, group=gr, p = 1, alpha = a)
    Resa=propd(counts=M, group=gr, p = 1, weighted = TRUE, alpha = a)

    #from these results, we need only theta and LRV:
    st=res@results[,c("lrv","theta")]
    stw=Res@results[,c("lrv","theta")]
    sta=resa@results[,c("lrv","theta")]
    stwa=Resa@results[,c("lrv","theta")]

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

    moda=dz*s2z/sta[,"lrv"]
    Fpmoda=(1-sta[,"theta"])*(dz+nce+nco)/((nce+nco)*sta[,"theta"]+moda)
    thetamoda=1/(1+Fpmoda)
    Fmoda=(nce+nco+dz-2)*Fpmoda

    #print("(4) weighted power-transformed moderated statistic")

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
      pd.nn@results$theta_mod,
      ref[[1]]
    )

    expect_equal(
      pd.nn@results$Fstat,
      ref[[5]]
    )

    expect_equal(
      pd.wn@results$theta_mod,
      ref[[2]]
    )

    expect_equal(
      pd.wn@results$Fstat,
      ref[[6]]
    )

    expect_equal(
      pd.na@results$theta_mod,
      ref[[3]]
    )

    expect_equal(
      pd.na@results$Fstat,
      ref[[7]]
    )

    expect_equal(
      pd.wa@results$theta_mod,
      ref[[4]]
    )

    expect_equal(
      pd.wa@results$Fstat,
      ref[[8]]
    )
  })
}
