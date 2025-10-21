library(testthat)
library(propr)
library(MASS)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("getCutoffFDR gets the correct cutoff", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get cutoff
    cutoff <- getCutoffFDR(pr, fdr = 0.05, window_size = 1)

    # check that the cutoff is correct
    expect_equal(round(cutoff, 4), 0.9738)  # for the moment the expected value is manually calculated
})

test_that("getCutoffFstat gets the correct cutoff when FDR correction is considered", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propd(x, as.character(y), p=10)
    pr <- updateF(pr, moderated=F)

    # get cutoff 
    cutoff_expected <- max(pr@results$theta[pr@results$FDR <= 0.05])
    cutoff_actual <- getCutoffFstat(pr, pval = 0.05, fdr = TRUE)

    # check that the cutoff is correct
    expect_equal(round(cutoff_actual, 6), round(cutoff_expected, 6))
})

test_that("getCutoffFstat gets the correct cutoff when FDR correction is not considered", {

    # get propr object and update cutoffs
    set.seed(0)
    pr <- propd(x, as.character(y), p=10)
    pr <- updateF(pr, moderated=F)

    # get cutoff 
    pval <- 0.05
    K <- length(unique(pr@group))
    N <- length(pr@group) + pr@dfz # population-level metric (i.e., N)
    Q <- stats::qf(pval, K - 1, N - K, lower.tail = FALSE)
    cutoff_expected <- (N - 2) / (Q + (N - 2))
    cutoff_actual <- getCutoffFstat(pr, fdr = FALSE)

    # check that the cutoff is correct
    expect_equal(round(cutoff_actual, 6), round(cutoff_expected, 6))
})