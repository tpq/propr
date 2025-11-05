library(testthat)
library(propr)
library(MASS)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("getAdjacencyFDR works properly for theta", {

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    adj_expected[propr:::getMatrix(pr) <= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(x)
    colnames(adj_expected) <- colnames(x)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjacencyFDR and getSignificantResultsFDR return coherent results", {

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjacencyFDR(pr)

    # get significant results
    results <- getSignificantResultsFDR(pr)

    # check that the values are correct
    for (i in 1:nrow(results)){
        expect_equal(adj[results[i,1], results[i,2]], 1)
    }
})

test_that("getAdjacencyFstat works properly", {

    for (fdr_adjusted in c(TRUE, FALSE)){

        # get propd object
        pr <- propd(x, as.character(y), p=10)
        pr <- updateF(pr, moderated=F)

        # get adjacency matrix
        adj <- getAdjacencyFstat(pr, fdr_adjusted=fdr_adjusted)

        # get expected adjacency matrix
        adj_expected <- matrix(0, nrow = ncol(x), ncol = ncol(x))
        adj_expected[propr:::getMatrix(pr) <= getCutoffFstat(pr, fdr_adjusted=fdr_adjusted)] <- 1
        adj_expected[diag(adj_expected)] <- 1
        rownames(adj_expected) <- colnames(x)
        colnames(adj_expected) <- colnames(x)

        # check that the values are correct
        expect_equal(adj, adj_expected)
    }
})

test_that("getAdjacencyFstat and getSignificantResultsFstat return coherent results", {

    for (fdr in c(TRUE, FALSE)){

        # get propd object
        pr <- propd(x, as.character(y), p=10)
        pr <- updateF(pr, moderated=F)

        # get adjacency matrix
        adj <- getAdjacencyFstat(pr, fdr=fdr)

        # get significant results
        results <- getSignificantResultsFstat(pr, fdr=fdr)

        # check that the values are correct
        for (i in 1:nrow(results)){
            expect_equal(adj[results[i,1], results[i,2]], 1)
        }
    }
})