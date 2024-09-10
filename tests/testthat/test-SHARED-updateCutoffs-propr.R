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

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)

  # get cutoffs
  cutoffs <- as.numeric( quantile(pr@matrix[lower.tri(pr@matrix)], probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(pr@fdr$cutoff, cutoffs)
})

test_that("updateCutoffs.propr properly calculates truecounts", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)

  # get truecounts
  truecounts1 <- sapply(
    pr@fdr$cutoff[pr@fdr$cutoff < 0], 
    function(cut) sum(pr@results$propr < cut)
  )
  truecounts2 <- sapply(
    pr@fdr$cutoff[pr@fdr$cutoff >= 0], 
    function(cut) sum(pr@results$propr > cut)
  )
  truecounts <- c(truecounts1, truecounts2)
  truecounts_manual <- c(0,1,2, 6,5,4,3,3,1,0)

  # check that truecounts are properly defined
  expect_equal(pr@fdr$truecounts, truecounts)
  expect_equal(pr@fdr$truecounts, truecounts_manual)
})

test_that("updateCutoffs.propr properly calculates randcounts", {

  # get propr object and update cutoffs
  set.seed(0)
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
      cutoff <- pr@fdr$cutoff[cut]
      if (cutoff >= 0) {
        randcounts[cut] <- randcounts[cut] + sum(pkt > cutoff)
      } else {
        randcounts[cut] <- randcounts[cut] + sum(pkt < cutoff)
      }
    }
  }
  randcounts <- randcounts / 10
  randcounts_manual <- c(0,0,0, 9.7,9.7,9.2,8.2,0,0,0)

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
