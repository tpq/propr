library(testthat)
library(propr)


test_that("propr - test that permute stay the same only when seed is given",{

  # define data matrix
  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)

  # test when seed is given
  set.seed(0)
  pr1 <- propr(X, metric = "pcor.bshrink", p=10)
  set.seed(0)
  pr2 <- propr(X, metric = "pcor.bshrink", p=10)
  expect_equal(
    pr1@permutes,
    pr2@permutes
  )

  # test when seed is not given
  pr1 <- propr(X, metric = "pcor.bshrink", p=10)
  pr2 <- propr(X, metric = "pcor.bshrink", p=10)
  expect_false(
    isTRUE(all.equal(
        pr1@permutes,
        pr2@permutes
    ))
  )
})

test_that("propd - test that permute stay the same only when seed is given",{

  # define data
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector
  v <- vegan::rda(log(x[,1]/x[,2]) ~ y)

  # test when seed is given
  set.seed(0)
  pd1 <- propd(x, as.character(y), p=10)
  set.seed(0)
  pd2 <- propd(x, as.character(y), p=10)
  expect_equal(
    pd1@permutes,
    pd2@permutes
  )

  # test when seed is not given
  pd1 <- propd(x, as.character(y), p=10)
  pd2 <- propd(x, as.character(y), p=10)
  expect_false(
    isTRUE(all.equal(
        pd1@permutes,
        pd2@permutes
    ))
  )
})

test_that("propr - test that permute conserves gene-wise or sample-wise totals", {
  # define data matrix
  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)

  # test feature-wise permutation
  set.seed(0)
  pr1 <- propr(X, metric = "pcor.bshrink", p=10, permutation_option = "feature-wise")
  expect_equal(
    as.vector(colSums(pr1@permutes[[1]])),
    as.vector(colSums(X))
  )

  # test sample-wise permutation
  set.seed(0)
  pr2 <- propr(X, metric = "pcor.bshrink", p=10, permutation_option = "sample-wise")
  expect_equal(
    as.vector(rowSums(pr2@permutes[[1]])),
    as.vector(rowSums(X))
  )
})