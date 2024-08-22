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

test_that("test that cutoffs are properly set  - for propr object", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)

  # get cutoffs
  cutoffs <- as.numeric( quantile(pr@matrix[lower.tri(pr@matrix)], probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(pr@fdr$cutoff, cutoffs)
})

test_that("test that the permutation tests work properly - for propr object", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)

  # get permuted values
  randcounts <- rep(0, 10)
  for (k in 1:10){
    ct.k <- pr@permutes[[k]]
    pr.k <- suppressMessages(propr(
      ct.k,
      pr@metric,
      ivar = pr@ivar,
      alpha = pr@alpha,
      p = 0
    ))
    pkt <- pr.k@results$propr
    for (cut in 1:length(pr@fdr$cutoff)){
        randcounts[cut] <- randcounts[cut] + propr:::countValuesBeyondThreshold(pkt, pr@fdr$cutoff[cut], pr@direct)
    }
  }
  randcounts <- randcounts / 10

  # check that the permutation tests work properly
  expect_equal(pr@fdr$randcounts, randcounts)
})

test_that("test that updateCutoffs can be reproducible when seed is set - for propr object", {
  
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
