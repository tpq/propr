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


test_that("getAdjacencyFDR returns the expected values for pcor.bshrink - clr", {

    # get propr object
    pr <- propr(X, metric = "pcor.bshrink", ivar='clr', p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix >= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR returns the expected values for pcor.bshrink - alr", {

    # get propr object
    pr <- propr(X, metric = "pcor.bshrink", ivar='alr', p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix >= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR returns the expected values for rho - clr", {

    # get propr object
    pr <- propr(X, metric = "rho", ivar='clr', p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix >= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR returns the expected values for rho - 5", {

    # get propr object
    pr <- propr(X, metric = "rho", ivar=5, p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix >= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR returns the expected values for phs - clr", {

    # get propr object
    pr <- propr(X, metric = "phs", ivar='clr', p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix <= getCutoffFDR(pr,)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR returns the expected values for phs - 5", {

    # get propr object
    pr <- propr(X, metric = "phs", ivar=5, p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=100)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    adj_expected[pr@matrix <= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(X)
    colnames(adj_expected) <- colnames(X)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR and getSignificantResultsFDR return coherent results", {
    
    for (metric in c('rho', 'phi', 'phs', 'pcor', 'pcor.bshrink')) { # pcor.shrink does not provide positive values for this dataset, and it gives error when tails = 'right'
        print(metric)

        # get propr object
        pr <- propr(X, metric=metric, p=10)

        # update FDR values
        pr <- updateCutoffs(pr, number_of_cutoffs=100)

        # get adjacency matrix
        adj <- getAdjacencyFDR(pr)

        # get significant results
        results <- getSignificantResultsFDR(pr)

        # check that the values are correct
        for (i in 1:nrow(results)){
            expect_equal(adj[results[i,1], results[i,2]], 1)
        }
    }
})