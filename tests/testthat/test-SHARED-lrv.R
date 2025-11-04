library(testthat)
library(propr)

# use iris data
counts <- iris[, 1:4]
counts <- as.matrix(counts)
group <- iris[, 5]

test_that("lrv_with_shrinkage and lrv match, when shrinkage is false", {

    # calculate lrv without shrinkage
    lrv1 <- propr:::lrv(counts, counts, FALSE, NA, counts, counts)
    lrv2 <- propr:::lrv_with_shrinkage(counts, shrink=FALSE)

    # check if they are equal
    expect_equal(lrv1, lrv2)
})


test_that("lrv using gpu and lrv match", {
    # calculate lrv 
    lrv1 <- propr:::lrv(counts, counts, FALSE, NA, counts, counts, use_gpu=TRUE)
    lrv2 <- propr:::lrv(counts, counts, FALSE, NA, counts, counts, use_gpu=FALSE)
    # check if they are equal
    # skip()
    expect_equal(lrv1, lrv2, tolerance = 1e-5)
})

test_that("lrv with_shrinkage and lrv do not match, when shrinkage is true", {
    # calculate lrv with shrinkage
    lrv1 <- propr:::lrv(counts, counts, FALSE, NA, counts, counts)
    lrv2 <- propr:::lrv_with_shrinkage(counts, shrink=TRUE)

    # check if they are not equal
    expect_false(identical(lrv1, lrv2))
})

test_that("lrv_with_shrinkage and lrv match for within group", {
    # define inputs
    ct <- as.matrix(counts)
    w <- ct
    groups <- lapply(unique(group), function(g) g == group)
    ngrp <- length(unique(group))

    # calculate lrv for group 1
    lrv1 <- propr:::lrv(ct[groups[[1]],], w[groups[[1]],], FALSE, NA, ct, w)
    lrv2 <- propr:::lrv_with_shrinkage(ct[groups[[1]],], shrink=FALSE)
    expect_equal(lrv1, lrv2)

    # calculate lrv for group 2
    lrv1 <- propr:::lrv(ct[groups[[2]],], w[groups[[2]],], FALSE, NA, ct, w)
    lrv2 <- propr:::lrv_with_shrinkage(ct[groups[[2]],], shrink=FALSE)
    expect_equal(lrv1, lrv2)
})