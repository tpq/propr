library(testthat)
library(propr)
library(corpcor)

# define data
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)
for (i in 1:4){
  X[sample(1:N, 10),i] <- 0
}

test_that("pcor.bshrink is correct when ivar is clr or alr",{

  for (ivar in c("clr", "alr")){

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    out <- propr:::pcor.bshrink(ct, ivar)
    mat <- out$matrix
    lambda <- out$lambda

    # compute pcor
    pr <- propr(X, metric = "pcor.bshrink", ivar=ivar)

    # expect counts to have zeros replaced
    expect_true(
      all(pr@counts == ct)
    )

    # NOTE that the data is not logratio transformed while computing pcor.bshrink
    # it is internally handled by covariance conversion.
    # so pr@logratio should be equal to pr@counts
    expect_true(
      all(pr@logratio == ct)
    )
  
    # expect computed coefficients are equal
    expect_true(
      all(round(pr@matrix, 8) == round(mat, 8))
    )

    # expect same lambda
    expect_equal(
      pr@lambda,
      lambda
    )

    # check dimensions are correct
    expect_equal(ncol(pr@matrix), 5)
    expect_equal(nrow(pr@matrix), 5)
  }
  
})

test_that("test that pcor.bshrink gives error when ivar is NA", {
    expect_error(
        propr(X, metric = "pcor.bshrink", ivar=NA)
    )
})

test_that("pcor.bshrink with alr and clr are the same", {

    # compute pcor with alr
    pr_alr <- suppressWarnings(propr(X, metric = "pcor.bshrink", ivar='alr'))
    pcor_alr <- getMatrix(pr_alr)[1:4, 1:4]

    # compute pcor with clr
    pr_clr <- propr(X, metric = "pcor.bshrink", ivar="clr")
    pcor_clr <- getMatrix(pr_clr)[1:4, 1:4]

    # expect that the coefficients are the same
    expect_equal(
        round(pcor_alr, 8),
        round(pcor_clr, 8)
    )

})
