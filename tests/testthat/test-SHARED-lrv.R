library(testthat)
library(propr)

# use iris data
counts <- iris[, 1:4]
counts <- as.matrix(counts)

test_that("lrv_with_shrinkage and lrv match, when shrinkage is false", {

    # calculate lrv without shrinkage
    lrv1 <- propr:::lrv(counts, counts, FALSE, NA, counts, counts)
    lrv2 <- propr:::lrv_with_shrinkage(counts, shrink=FALSE)

    # check if they are equal
    expect_equal(lrv1, lrv2)
})

test_that("lrv with_shrinkage and lrv do not match, when shrinkage is true", {
    # calculate lrv with shrinkage
    lrv1 <- propr:::lrv(counts, counts, TRUE, NA, counts, counts)
    lrv2 <- propr:::lrv_with_shrinkage(counts, shrink=TRUE)

    # check if they are not equal
    expect_false(identical(lrv1, lrv2))
})
