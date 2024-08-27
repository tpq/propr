library(testthat)
library(propr)

message_test <- function(title) {
    message(
        "==========================================================\n", 
        "....Running test: ", title, "\n")
}

# define data matrix
set.seed(123)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)


test_that("getAdjFDR returns the expected values for pcor.bshrink", {

    message_test("getAdjFDR returns the expected values for pcor.bshrink")

    # get propr object
    pr <- propr(X, metric = "pcor.bshrink", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr, consider_negative_values=TRUE)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix <= getCutoffFDR(pr, positive = FALSE) | pr@matrix >= getCutoffFDR(pr, positive = TRUE)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFDR returns the expected values for rho", {

    message_test("getAdjFDR returns the expected values for rho")

    # get propr object
    pr <- propr(X, metric = "rho", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr, consider_negative_values=FALSE)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix >= getCutoffFDR(pr, positive = TRUE)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFDR returns the expected values for phs", {

    message_test("getAdjFDR returns the expected values for phs")

    # get propr object
    pr <- propr(X, metric = "phs", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr, consider_negative_values=FALSE)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix <= getCutoffFDR(pr, positive = TRUE)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFDR and getSignificantResultsFDR return coherent results", {
    
    message_test("getAdjFDR and getSignificantResultsFDR return coherent results")

    # get propr object
    set.seed(0)
    pr <- propr(X, metric = "pcor.bshrink", p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr, consider_negative_values=TRUE)

    # get significant results
    results <- getSignificantResultsFDR(pr, consider_negative_values=TRUE)

    # check that the values are correct
    for (i in 1:nrow(results)){
        expect_equal(adj[results[i,1], results[i,2]], 1)
    }
})