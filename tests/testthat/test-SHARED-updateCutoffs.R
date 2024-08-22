library(testthat)
library(propr)
library(MASS)

test_that("test that cutoffs are properly set  - for propr object", {

  # define data matrix
  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)

  # get cutoffs
  cutoffs <- as.numeric( quantile(pr@matrix[lower.tri(pr@matrix)], probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(pr@fdr$cutoff, cutoffs)
})

test_that("test that updateCutoffs can be reproducible when seed is set - for propr object", {
  
    # define data matrix
    N <- 100
    a <- seq(from = 5, to = 15, length.out = N)
    b <- a * rnorm(N, mean = 1, sd = 0.1)
    c <- rnorm(N, mean = 10)
    d <- rnorm(N, mean = 10)
    e <- rep(10, N)
    X <- data.frame(a, b, c, d, e)
  
    # get propr object and update cutoffs
    set.seed(0)
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    set.seed(0)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
  
    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})

test_that("test that updateCutoffs will give different permutation results when seed is not set - for propr object", {
    
    # define data matrix
    N <- 100
    a <- seq(from = 5, to = 15, length.out = N)
    b <- a * rnorm(N, mean = 1, sd = 0.1)
    c <- rnorm(N, mean = 10)
    d <- rnorm(N, mean = 10)
    e <- rep(10, N)
    X <- data.frame(a, b, c, d, e)
    
    # get propr object and update cutoffs
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
    
    # check that fdr are different
    expect_false(isTRUE(all.equal(pr1@fdr, pr2@fdr)))

    # check that at least the cutoffs are the same
    expect_equal(pr1@fdr$cutoff, pr2@fdr$cutoff)
})

test_that("test that updateCutoffs works when ncores > 1 - for propr object", {

    # define data matrix
    N <- 100
    a <- seq(from = 5, to = 15, length.out = N)
    b <- a * rnorm(N, mean = 1, sd = 0.1)
    c <- rnorm(N, mean = 10)
    d <- rnorm(N, mean = 10)
    e <- rep(10, N)
    X <- data.frame(a, b, c, d, e)

    # get propr object and update cutoffs
    set.seed(0)
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10, ncores=1)
    set.seed(0)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10, ncores=2)

    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})

test_that("test that cutoffs are properly set - for propd object", {

  # define data
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector

  # get propd object and update cutoffs
  pd <- propd(x, as.character(y), p=10)
  pd <- updateCutoffs(pd, number_of_cutoffs=10)

  # get cutoffs
  cutoffs <- as.numeric( quantile(pd@results$theta, probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(pd@fdr$cutoff, cutoffs)
})

test_that("test that updateCutoffs can be reproducible when seed is set - for propr object", {
  
    # define data
    data(crabs)
    x <- crabs[,4:8]  # data matrix with 5 variables
    y <- crabs[,1]    # group vector
    
    # get propd object and update cutoffs
    set.seed(0)
    pr1 <- propd(x, as.character(y), p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    set.seed(0)
    pr2 <- propd(x, as.character(y), p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
  
    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})

test_that("test that updateCutoffs will give different permutation results when seed is not set - for propd object", {
    
    # define data
    data(crabs)
    x <- crabs[,4:8]  # data matrix with 5 variables
    y <- crabs[,1]    # group vector
    
    # get propd object and update cutoffs
    pr1 <- propd(x, as.character(y), p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    pr2 <- propd(x, as.character(y), p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
    
    # check that fdr are different
    expect_false(isTRUE(all.equal(pr1@fdr, pr2@fdr)))

    # check that at least the cutoffs are the same
    expect_equal(pr1@fdr$cutoff, pr2@fdr$cutoff)
})

test_that("test that updateCutoffs works when ncores > 1 - for propd object", {
  
    # define data
    data(crabs)
    x <- crabs[,4:8]  # data matrix with 5 variables
    y <- crabs[,1]    # group vector

    # get propd object and update cutoffs
    set.seed(0)
    pr1 <- propd(x, as.character(y), p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10, ncores=1)
    set.seed(0)
    pr2 <- propd(x, as.character(y), p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10, ncores=2)

    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})
