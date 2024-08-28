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

test_that("pcor.shrink is the same when ivar is 5, with or without the reference column", {

  # compute lr
  ct <- simple_zero_replacement(X)
  lr <- logratio_without_alpha(ct, 5)

  # compute pcor.shrink
  pr1 <- suppressWarnings(propr(lr, metric = "pcor.shrink", ivar=NA))

  # compute pcor.shrink without 5th column
  pr2 <- suppressWarnings(propr(lr[1:4], metric = "pcor.shrink", ivar=NA))

  # expect computed coefficients are equal
  expect_true(
    all(round(pr1@matrix, 8)[1:4,1:4] == round(pr2@matrix, 8)[1:4,1:4])
  )

  # check shrinkage is the same
  expect_equal(
    pr1@lambda,
    pr2@lambda
  )

})

test_that("pcor.shrink is correct when ivar is given", {

  d1 <- list("clr", 5, c(1,3))
  d2 <- list(c(1:5), 5, c(1,3))

  for (i in 1:length(d1)){

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    lr <- logratio_without_alpha(ct, d2[[i]])
    cov <- suppressWarnings(cov.shrink(lr))
    mat <- cor2pcor(cov)
    mat <- matrix(mat, ncol=ncol(lr), nrow=ncol(lr))
    class(mat) <- "matrix"
    lambda <- attr(cov, "lambda")

    # compute pcor
    pr <- suppressWarnings(propr(X, metric = "pcor.shrink", ivar=d1[[i]]))

    # expect counts to have zeros replaced
    expect_true(
      all(pr@counts == ct)
    )

    # expect that logratio are equal
    expect_true(
      all(pr@logratio == lr)
    )

    # expect computed coefficients are equal
    expect_true(
      all(round(pr@matrix, 8) == round(mat, 8), na.rm=T)
    )

    # check shrinkage is not applied
    expect_equal(
      pr@lambda,
      lambda
    )
  }
})


test_that("pcor.shrink is correct when ivar is NA using previously transformed data", {

  d1 <- list("clr", 5, c(1,3))
  d2 <- list(c(1:5), 5, c(1,3))

  for (i in 1:length(d1)){

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    lr <- logratio_without_alpha(ct, d2[[i]])
    cov <- suppressWarnings(cov.shrink(lr))
    mat <- cor2pcor(cov)
    mat <- matrix(mat, ncol=ncol(lr), nrow=ncol(lr))
    class(mat) <- "matrix"
    lambda <- attr(cov, "lambda")

    # compute pcor
    pr <- suppressWarnings(propr(lr, metric = "pcor.shrink", ivar=NA))

    # expect counts to contain the previously transformed data
    expect_true(
      all(pr@counts == lr)
    )

    # expect that the logratio slot also contains the exact input data
    expect_true(
      all(pr@logratio == lr)
    )

    # expect computed coefficients are equal
    expect_true(
      all(round(pr@matrix, 8) == round(mat, 8), na.rm=T)
    )

    # check shrinkage is not applied
    expect_equal(
      pr@lambda,
      lambda
    )
  }
})

test_that("pcor.shrink with alr and clr are not the same", {

    # compute pcor.shrink with alr
    pr <- suppressWarnings(propr(X, metric = "pcor.shrink", ivar=5))
    pcor_alr <- getMatrix(pr)[1:4, 1:4]

    # compute pcor.shrink with clr
    pr <- propr(X, metric = "pcor.shrink", ivar="clr")
    pcor_clr <- getMatrix(pr)[1:4, 1:4]

    expect_false(
      all(round(pcor_alr, 8) == round(pcor_clr, 8))
    )

})
