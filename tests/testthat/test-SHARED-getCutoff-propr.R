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

test_that("test that getCutoff gets the correct cutoff - for positive propr values", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propr(X, metric = "pcor.bshrink", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1, positive = TRUE)

    # check that the cutoff is correct
    expect_equal(round(cutoff, 4), 0.6573)  # for the moment the expected value is manually calculated
})

test_that("test that getCutoff gets the correct cutoff - for negative propr values", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propr(X, metric = "pcor.bshrink", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1, positive = FALSE)

    # check that the cutoff is correct
    expect_equal(round(cutoff, 4), -0.0313)  # for the moment the expected value is manually calculated
})