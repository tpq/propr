library(testthat)
library(propr)
library(MASS)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("test that getCutoff gets the correct cutoff", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)

    # check that the cutoff is correct
    expect_equal(round(cutoff, 4), 0.9738)  # for the moment the expected value is manually calculated
})