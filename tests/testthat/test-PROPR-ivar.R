library(testthat)
library(propr)

data(mtcars)

test_that("the data is properly handled when ivar is NA", {

    # compute propr object
    pr <- propr(mtcars, ivar=NA)

    # check the counts remain the same as the original input
    expect_equal(
        pr@counts, 
        mtcars
    )

    # check the logratios are the same
    expect_equal(
        pr@logratio, 
        mtcars
    )
})

test_that("the data is properly handled when ivar is clr", {

    # compute propr object
    pr <- propr(mtcars, ivar="clr")

    # compute expected data
    ct <- simple_zero_replacement(mtcars)
    clr <- logratio_without_alpha(ct, c(1:ncol(ct)))

    # check the counts remain the same as the zero-replaced input
    expect_equal(
        pr@counts, 
        ct
    )

    # check the logratios are the same
    expect_equal(
        pr@logratio, 
        clr
    )
})

test_that("the data is properly handled when ivar is 3", {

    # compute propr object
    pr <- propr(mtcars, ivar=3)

    # compute expected data
    ct <- simple_zero_replacement(mtcars)
    alr <- logratio_without_alpha(ct, 3)

    # check the counts remain the same as the original input
    expect_equal(
        pr@counts, 
        ct
    )

    # check the logratios are the same
    expect_equal(
        pr@logratio, 
        alr
    )
})

test_that("the data is properly handled when ivar is 1,6", {

    # compute propr object
    pr <- propr(mtcars, ivar=c(1,6))

    # compute expected data
    ct <- simple_zero_replacement(mtcars)
    alr <- logratio_without_alpha(ct, c(1,6))

    # check the counts remain the same as the original input
    expect_equal(
        pr@counts, 
        ct
    )

    # check the logratios are the same
    expect_equal(
        pr@logratio, 
        alr
    )
})

test_that("pearson correlation is correct when ivar is NA", {

    # compute propr object
    cor_propr <- propr(mtcars, metric = "cor", ivar=NA)@matrix
    
    # get correlation using cor
    cor_cor <- cor(mtcars, method = "pearson")

    expect_equal(
        round(cor_propr, 6), 
        round(cor_cor, 6)
    )
})

test_that("pearson correlation is correct when ivar is clr", {

    # get correlation using propr
    pr <- propr(mtcars, metric = "cor", ivar="clr")

    # get correlation using cor
    ct <- simple_zero_replacement(mtcars)
    clr <- logratio_without_alpha(ct, c(1:ncol(ct)))
    ccor <- cor(clr, method = "pearson")

    # check the correlations are the same
    expect_equal(
        round(pr@matrix, 6), 
        round(ccor, 6)
    )
})

test_that("pearson correlation is correct when ivar is 1", {

    # get correlation using propr
    pr <- suppressWarnings(propr(mtcars, metric = "cor", ivar=1))

    # get correlation using cor
    ct <- simple_zero_replacement(mtcars)
    alr <- logratio_without_alpha(ct, 1)
    ccor <- suppressWarnings(cor(alr, method = "pearson"))

    # check the correlations are the same
    expect_equal(
        round(pr@matrix, 6), 
        round(ccor, 6)
    )
})

test_that("pearson correlation is correct when ivar is 1,3", {

    # get correlation using propr
    pr <- suppressWarnings(propr(mtcars, metric = "cor", ivar=c(1,3)))

    # get correlation using cor
    ct <- simple_zero_replacement(mtcars)
    alr <- logratio_without_alpha(ct, c(1,3))
    ccor <- suppressWarnings(cor(alr, method = "pearson"))

    # check the correlations are the same
    expect_equal(
        round(pr@matrix, 6), 
        round(ccor, 6)
    )
})
