library(testthat)
library(propr)

data(mtcars)


test_that("pearson correlation is correct when ivar is NA", {

    # get correlation using propr
    cor_propr <- propr(mtcars, metric = "cor", ivar=NA)@matrix
    cor_propr <- round(cor_propr, 6)
    
    # get correlation using cor
    cor_cor <- cor(mtcars, method = "pearson")
    cor_cor <- round(cor_cor, 6)

    expect_equal(cor_propr, cor_cor)
})

test_that("pearson correlation is correct when ivar is clr", {

    # get correlation using propr
    pr <- propr(mtcars, metric = "cor", ivar="clr")

    # get correlation using cor
    ct <- simple_zero_replacement(mtcars)
    clr <- t(apply(ct, 1, function(x) log(x) - mean(log(x))))
    ccor <- cor(clr, method = "pearson")

    # check the counts are the same
    expect_equal(
        pr@counts, 
        ct
    )
    
    # check the logratios are the same
    expect_equal(
        as.matrix(round(pr@logratio, 6)), 
        as.matrix(round(clr, 6))
    )

    # check the correlations are the same
    expect_equal(
        round(pr@matrix, 6), 
        round(ccor, 6)
    )
})

test_that("pearson correlation is correct when ivar is the first variable", {

    ref = 1

    # get correlation using propr
    pr <- propr(mtcars, metric = "cor", ivar=ref)

    # get correlation using cor
    ct <- simple_zero_replacement(mtcars)
    alr <- log(ct) - log(ct[,ref])
    ccor <- cor(alr, method = "pearson")

    # check the counts are the same
    expect_equal(
        pr@counts, 
        ct
    )
    
    # check the logratios are the same
    expect_equal(
        as.matrix(round(pr@logratio, 6)), 
        as.matrix(round(alr, 6))
    )

    # check the correlations are the same
    expect_equal(
        round(pr@matrix, 6), 
        round(ccor, 6)
    )
})


