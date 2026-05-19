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

test_that("test that getCutoff gets the correct cutoff - sample wise", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propr(X, metric = "pcor.bshrink", p=10, permutation_option = "sample-wise")
    pr <- updateCutoffs(pr, number_of_cutoffs=10, tails='right')

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)

    # check that cutoff is correct
    expect_equal(round(cutoff, 4), 0.6573)  # for the moment the expected value is manually calculated
})

test_that("test that getCutoff gets the correct cutoff - feature wise", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propr(X, metric = "pcor.bshrink", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10, tails='right')

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)

    # check that cutoff is correct
    expect_equal(round(cutoff, 4), 0.7762)  # for the moment the expected value is manually calculated
})

test_that("test that moving average works as expected", {

    x <- c(NA, 1, 1, 2, NA, 2, 3, 3, 3, 4, 4, 5)

    # define expected moving average
    y1 <- x 
    y2 <- c(NA, 1, 1.5, 2, NA, 2.5, 3, 3, 3.5, 4, 4.5, 5)
    y3 <- c(NA, 1, 1.33, 1.5, NA, 2.5, 2.67, 3, 3.33, 3.67, 4.33, 4.5)
    y4 <- c(NA, 1.33, 1.33, 1.67, NA, 2.67, 2.75, 3.25, 3.5, 4, 4.33, 4.5)

    # check
    expect_equal(round(propr:::getMovingAverage(x, 1),2), round(y1,2))
    expect_equal(round(propr:::getMovingAverage(x, 2),2), round(y2,2))
    expect_equal(round(propr:::getMovingAverage(x, 3),2), round(y3,2))
    expect_equal(round(propr:::getMovingAverage(x, 4),2), round(y4,2))
})