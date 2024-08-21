library(testthat)
library(propr)


test_that("test that updateCutoff work for propr",{

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
  pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
  set.seed(0)
  pr2 <- propr(X, metric = "pcor.bshrink", p=10)
  pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
  expect_equal(
    pr1@fdr,
    pr2@fdr
  )

  # test when seed is not given
  pr1 <- propr(X, metric = "pcor.bshrink", p=10)
  pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
  pr2 <- propr(X, metric = "pcor.bshrink", p=10)
  pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
  expect_false(
    isTRUE(all.equal(
        pr1@fdr,
        pr2@fdr
    ))
  )

  # test that cutoffs are properly defined
  expect_equal(
    pr1@fdr$cutoff,
    pr2@fdr$cutoff
  )
  expect_equal(
    pr1@fdr$cutoff,
    as.vector( quantile(pr1@matrix[lower.tri(pr1@matrix)], probs = seq(0, 1, length.out = 10)) )
  )
})

test_that("test that updateCutoff work for propd",{

  # define data
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector
  v <- vegan::rda(log(x[,1]/x[,2]) ~ y)

  # test when seed is given
  set.seed(0)
  pd1 <- propd(x, as.character(y), p=10)
  pd1 <- updateCutoffs(pd1, number_of_cutoffs=10)
  set.seed(0)
  pd2 <- propd(x, as.character(y), p=10)
  pd2 <- updateCutoffs(pd2, number_of_cutoffs=10)
  expect_equal(
    pd1@fdr,
    pd2@fdr
  )

  # test when seed is not given
  # for the moment, both fdr are the same, 
  # because it is extremely difficult to get 
  # some high theta values by change with this iris data
  # pd1 <- propd(x, as.character(y), p=10)
  # pd1 <- updateCutoffs(pd1, number_of_cutoffs=10)
  # pd2 <- propd(x, as.character(y), p=10)
  # pd2 <- updateCutoffs(pd2, number_of_cutoffs=10)
  # expect_false(
  #   isTRUE(all.equal(
  #       pd1@fdr,
  #       pd2@fdr
  #   ))
  # )

  # check that cutoffs are properly defined
  expect_equal(
    pd1@fdr$cutoff,
    pd2@fdr$cutoff
  )
  expect_equal(
    pd1@fdr$cutoff,
    as.vector( quantile(pd1@results$theta, probs = seq(0, 1, length.out = 10)) )
  )
})

