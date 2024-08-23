library(testthat)
library(propr)
library(MASS)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("test that getResults works as expected", {
  
  # get propr object
  pr <- propd(x, as.character(y), p=10)
  
  # get results
  results <- getResults(pr)
  
  # check that the values are correct
  expect_equal(pr@results[,-c(1:2)], results[,-c(1:2)])

  # check that the variable names are corretly replaced
  expect_equal(results[,1], colnames(x)[pr@results[,1]])
  expect_equal(results[,2], colnames(x)[pr@results[,2]])

  # check that the order of pairs are as expected
  ord <- list(c(2,1), c(3,1), c(3,2), c(4,1), c(4,2), c(4,3), c(5,1), c(5,2), c(5,3), c(5,4))
  for (i in 1:10) {
    expect_equal(results[i,1], colnames(x)[ord[[i]][1]])
    expect_equal(results[i,2], colnames(x)[ord[[i]][2]])
  }
})

test_that("test that getSignificantResultsFDR works as expected", {
  
  # get propr object
  pr <- propd(x, as.character(y), p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)
  
  # get expected results
  cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)
  expected <- pr@results[which(pr@results$theta <= cutoff),]
  
  # get significant results
  results <- getSignificantResultsFDR(pr, fdr = 0.05)
  
  # check that the values are correct
  expect_equal(results$theta, expected$theta)
})