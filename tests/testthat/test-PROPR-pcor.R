library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

test_that("pcor is correct when ivar is clr",{
  
    # compute pcor
    pr <- propr(X, metric = "pcor", ivar="clr")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(pr@counts)
    clr <- t(apply(ct, 1, function(x) log(x) - mean(log(x))))
    cov <- cov(clr)
    mat <- corpcor::cor2pcor(cov)
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    expect_equal(
      pr@counts,
      ct
    )
  
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )
  
  
})

test_that("pcor is correct when ivar is NA", {

  # compute pcor
  pr <- propr(X, metric = "pcor", ivar=NA)

  # compute pcor manually
  cov <- cov(X)
  mat <- corpcor::cor2pcor(cov)
  class(mat) <- "matrix"
  colnames(mat)  <- colnames(X)
  rownames(mat)  <- colnames(X)

  expect_equal(
    round(pr@matrix, 8),
    round(mat, 8)
  )

})

test_that("pcor with alr and clr are the same", {

  # compute pcor with alr
  pr <- propr(X, metric = "pcor", ivar=5)
  pcor_alr <- getMatrix(pr)[1:4, 1:4]

  # compute pcor with clr
  pr <- propr(X, metric = "pcor", ivar="clr")
  pcor_clr <- getMatrix(pr)[1:4, 1:4]

  expect_equal(
    round(pcor_alr, 8),
    round(pcor_clr, 8)
  )

})

test_that("pcor.bshrink with alr and clr are the same", {

  # compute pcor.bshrink with alr
  pr <- propr(X, metric = "pcor.bshrink", ivar="alr")
  pcor_alr <- getMatrix(pr)[1:4, 1:4]

  # compute pcor.bshrink with clr
  pr <- propr(X, metric = "pcor.bshrink", ivar="clr")
  pcor_clr <- getMatrix(pr)[1:4, 1:4]

  expect_equal(
    round(pcor_alr, 8),
    round(pcor_clr, 8)
  )

})
