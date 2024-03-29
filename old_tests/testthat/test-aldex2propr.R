if(requireNamespace("ALDEx2", quietly = TRUE) &
   requireNamespace("Biobase", quietly = TRUE)){

  if(as.numeric(gsub("\\.", "", packageVersion("ALDEx2"))) >= 1101){

    library(propr)
    library(ALDEx2)

    data(mail)
    x <- aldex.clr(as.data.frame(t(mail)), conds = rep("A", 5))
    y <- aldex.clr(as.data.frame(t(mail)), conds = rep("A", 5), denom = 2)

    data(marg.abs)
    z <- ALDEx2:::iqlr.features(t(marg.abs[, 1:50]), conds = rep("A", 16))

    propr.phisym <- function (X){

      Cov    <- stats::var(X)
      tmp    <- 2 * Cov / outer(diag(Cov), diag(Cov), "+")
      return((1-tmp)/(1+tmp))
    }

    codaSeq.phi <- function(aldex.clr){

      sym.phi <- propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
                                       function(y){y[,1]})))

      for(i in 2:numMCInstances(aldex.clr)){
        sym.phi <- sym.phi + propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
                                                   function(y){y[,i]})))
      }

      lt <- which(col(sym.phi)<row(sym.phi), arr.ind=FALSE)
      lt.ind <- which(col(sym.phi)<row(sym.phi), arr.ind=TRUE)
      sma.df <- data.frame(row=factor(rownames(sym.phi)[lt.ind[,"row"]]),
                           col=factor(colnames(sym.phi)[lt.ind[,"col"]]), stringsAsFactors=FALSE)
      sma.df$phi <- sym.phi[lt] /  numMCInstances(aldex.clr)
      return(sma.df)
    }

    test_that("aldex2propr matches codaSeq implementation", {

      expect_equal(
        as.vector(propr.phisym(propr:::clrRcpp(mail[]))),
        as.vector(phis(mail)@matrix)
      )

      expect_equal(
        sort(as.vector(codaSeq.phi(x)$phi)),
        sort(propr:::lltRcpp(aldex2propr(x, how = "phis")@matrix))
      )
    })

    test_that("aldex2propr handles ivar correctly", {

      expect_equal(
        aldex2propr(y, how = "phi")@matrix[, 2],
        c(Inf, 0, Inf, Inf)
      )

      expect_equal(
        aldex2propr(y, how = "rho")@matrix[, 2],
        c(0, 1, 0, 0)
      )

      expect_equal(
        aldex2propr(y, how = "phs")@matrix[, 2],
        c(Inf, 0, Inf, Inf)
      )
    })

    test_that("iqlr matches ALDEx2 implementation", {

      expect_equal(
        perb(marg.abs[, 1:50] + .5, ivar = "iqlr")@matrix,
        perb(marg.abs[, 1:50] + .5, ivar = z)@matrix
      )
    })
  }
}
