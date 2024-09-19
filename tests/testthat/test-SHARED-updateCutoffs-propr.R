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

test_that("updateCutoffs.propr properly set up cutoffs", {

  # get propr object
  pr <- propr(X, metric = "pcor.bshrink", p=10)

  # get cutoffs
  values <- pr@matrix[lower.tri(pr@matrix)]
  cutoffs_right <- as.numeric( quantile(values[values >= 0], probs = seq(0, 1, length.out = 10)) )
  cutoffs_both <- as.numeric( quantile(abs(values), probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10, tails="right")@fdr$cutoff,
    cutoffs_right
  )
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10, tails="both")@fdr$cutoff,
    cutoffs_both
  )
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10)@fdr$cutoff,
    cutoffs_both
  )
})

test_that("updateCutoffs.propr properly calculates truecounts", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr_right <- updateCutoffs(pr, number_of_cutoffs=10, tails='right')
  pr_both <- updateCutoffs(pr, number_of_cutoffs=10, tails='both')
  pr_default <- updateCutoffs(pr, number_of_cutoffs=10)

  # get truecounts
  truecounts_right <- sapply(
    pr@fdr$cutoff, 
    function(cut) sum(pr@results$propr > cut)
  )
  truecounts_both <- sapply(
    pr@fdr$cutoff, 
    function(cut) sum(abs(pr@results$propr) > cut)
  )

  # check that truecounts are properly defined
  expect_true(all(pr_right@fdr$truecounts == truecounts_right))
  expect_true(all(pr_both@fdr$truecounts == truecounts_both))
  expect_true(all(pr_default@fdr$truecounts == truecounts_both))
})

test_that("updateCutoffs.propr properly calculates randcounts", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10, tails='both')

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
    pkt <- abs(pkt)
    for (cut in 1:length(pr@fdr$cutoff)){
      cutoff <- pr@fdr$cutoff[cut]
      randcounts[cut] <- randcounts[cut] + sum(pkt > cutoff)
    }
  }
  randcounts <- randcounts / 10

  # check that the permutation tests work properly
  expect_equal(pr@fdr$randcounts, randcounts)
})

test_that("updateCutoffs.propr is reproducible when seed is set", {
  
    # get propr object and update cutoffs
    set.seed(0)
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    set.seed(0)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
    pr3 <- propr(X, metric = "pcor.bshrink", p=10)
    pr3 <- updateCutoffs(pr3, number_of_cutoffs=10)
  
    # check that fdr are the same only when seed is the same
    expect_equal(pr1@fdr, pr2@fdr)
    expect_false(isTRUE(all.equal(pr1@fdr, pr3@fdr)))
})

test_that("updateCutoffs.propr works when ncores > 1", {

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
