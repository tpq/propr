library(testthat)
library(propr)

# define data matrix
set.seed(123)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

test_that("test that getResults works as expected", {
  
  # get propr object
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  
  # get results
  results <- getResults(pr)
  
  # check that the values are correct
  expect_equal(pr@results[,c(3:7)], results[,c(3:7)])

  # check that the variable names are corretly replaced
  expect_equal(results[,1], colnames(X)[pr@results[,1]])
  expect_equal(results[,2], colnames(X)[pr@results[,2]])

  # check that the order of pairs are as expected
  ord <- list(c(2,1), c(3,1), c(3,2), c(4,1), c(4,2), c(4,3), c(5,1), c(5,2), c(5,3), c(5,4))
  for (i in 1:10) {
    expect_equal(results[i,1], colnames(X)[ord[[i]][1]])
    expect_equal(results[i,2], colnames(X)[ord[[i]][2]])
  }
})

test_that("test that getSignificantResultsFDR works as expected - pcor.bshrink", {
  
  # get propr object
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)
  
  # get expected results
  cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)
  expected <- pr@results$propr[which(abs(pr@results$propr) >= cutoff)]

  # get significant results
  results <- getSignificantResultsFDR(pr, fdr = 0.05)

  # check that the values are correct
  expect_equal(results$propr, expected)
})

test_that("test that getSignificantResultsFDR works as expected - rho", {
  
  # get propr object
  pr <- propr(X, metric = "rho", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)
  
  # get expected results
  cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)
  expected <- pr@results$propr[which(pr@results$propr >= cutoff)]

  # get significant results
  results <- getSignificantResultsFDR(pr, fdr = 0.05)

  # check that the values are correct
  expect_equal(results$propr, expected)
})