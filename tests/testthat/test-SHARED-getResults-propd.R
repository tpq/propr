library(testthat)
library(propr)
library(MASS)

message_test <- function(title) {
    message(
        "==========================================================\n", 
        "....Running test: ", title, "\n")
}

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("getResults works as expected", {

  message_test("getResults works as expected")
  
  # get propd object
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

test_that("getSignificantResultsFDR works as expected", {

  message_test("getSignificantResultsFDR works as expected")
  
  # get propd object
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

test_that("getSignificantResultsFstat works as expected when using theoretical Fstat cutoff", {

  message_test("getSignificantResultsFstat works as expected when using theoretical Fstat cutoff")

  # get propd object
  pr <- propd(x, as.character(y), p=10)
  pr <- updateF(pr)

  # get Fstat cutoff
  pval <- 0.05
  K <- length(unique(pr@group))
  N <- length(pr@group) + pr@dfz # population-level metric (i.e., N)
  Q <- stats::qf(pval, K - 1, N - K, lower.tail = FALSE)
  cutoff <- (N - 2) / (Q + (N - 2))

  # expect that the Fstat values are smaller or equal than the cutoff
  expect_true(all(getSignificantResultsFstat(pr, fdr=F)$theta <= cutoff))

})

test_that("getSignificantResultsFstat works as expected when using FDR corrected values", {

  message_test("getSignificantResultsFstat works as expected when using FDR corrected values")

  # get propd object
  pr <- propd(x, as.character(y), p=10)
  pr <- updateF(pr)

  expect_true(all(getSignificantResultsFstat(pr, fdr=T)$fdr <= 0.05))
  
})