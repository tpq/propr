library(testthat)
library(propr)
library(MASS)

message_test <- function(title) {
    message(
        "==========================================================\n", 
        "....Running test: ", title, "\n")
}

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("getAdjFDR works properly for theta", {

    message_test("getAdjFDR works properly for theta")

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    adj_expected[propr:::get_theta_matrix(pr) <= getCutoffFDR(pr)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(x)
    colnames(adj_expected) <- colnames(x)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFDR and getSignificantResultsFDR return coherent results", {

    message_test("getAdjFDR and getSignificantResultsFDR return coherent results")

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateCutoffs(pr, number_of_cutoffs=10)

    # get adjacency matrix
    adj <- getAdjFDR(pr)

    # get significant results
    results <- getSignificantResultsFDR(pr)

    # check that the values are correct
    for (i in 1:nrow(results)){
        expect_equal(adj[results[i,1], results[i,2]], 1)
    }
})

test_that("getAdjFstat works properly when considering theoretical F-stat cutoff", {

    message_test("getAdjFstat works properly when considering theoretical F-stat cutoff")

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateF(pr)

    # get adjacency matrix
    adj <- getAdjFstat(pr, fdr=F)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    adj_expected[propr:::get_theta_matrix(pr) <= getCutoffFstat(pr, fdr=F)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(x)
    colnames(adj_expected) <- colnames(x)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFstat works properly when considering FDR corrected Fstatistics", {

    message_test("getAdjFstat works properly when considering FDR corrected Fstatistics")

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateF(pr)

    # get adjacency matrix
    adj <- getAdjFstat(pr, fdr=T)

    # get expected adjacency matrix
    adj_expected <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    adj_expected[propr:::get_theta_matrix(pr) <= getCutoffFstat(pr, fdr=T)] <- 1
    adj_expected[diag(adj_expected)] <- 1
    rownames(adj_expected) <- colnames(x)
    colnames(adj_expected) <- colnames(x)

    # check that the values are correct
    expect_equal(adj, adj_expected)
})

test_that("getAdjFstat and getSignificantResultsFstat return coherent results when considering theoretical F-stat cutoff", {

    message_test("getAdjFstat and getSignificantResultsFstat return coherent results when considering theoretical F-stat cutoff")

    # get propd object
    pr <- propd(x, as.character(y), p=10)
    pr <- updateF(pr)

    # get adjacency matrix
    adj <- getAdjFstat(pr, fdr=F)

    # get significant results
    results <- getSignificantResultsFstat(pr, fdr=F)

    # check that the values are correct
    for (i in 1:nrow(results)){
        expect_equal(adj[results[i,1], results[i,2]], 1)
    }
})

# test_that("getAdjFstat and getSignificantResultsFstat return coherent results when considering FDR corrected Fstatistics", {

#     message_test("getAdjFstat and getSignificantResultsFstat return coherent results when considering FDR corrected Fstatistics")

#     # get propd object
#     pr <- propd(x, as.character(y), p=10)
#     pr <- updateF(pr)

#     # get adjacency matrix
#     adj <- getAdjFstat(pr, fdr=T)

#     # get significant results
#     results <- getSignificantResultsFstat(pr, fdr=T)

#     # check that the values are correct
#     for (i in 1:nrow(results)){
#         expect_equal(adj[results[i,1], results[i,2]], 1)
#     }
# })