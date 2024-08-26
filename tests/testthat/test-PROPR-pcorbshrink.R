library(testthat)
library(propr)

message_test <- function(title) {
    message(
        "==========================================================\n", 
        "....Running test: ", title, "\n")
}

# define data
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

test_that("pcor.bshrink is correct when ivar is clr",{

    message_test("pcor.bshrink is correct when ivar is clr")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    out <- propr:::basis_shrinkage(ct)
    mat <- out$pcor
    lambda <- out$lambda
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- propr(X, metric = "pcor.bshrink", ivar="clr")

    # expect counts to have zeros replaced
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # NOTE that the data is not logratio transformed while computing pcor.bshrink
    # it is internally handled by covariance conversion.
    # so pr@logratio should be equal to pr@counts
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(pr@counts)
    )
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expect same lambda
    expect_equal(
      pr@lambda,
      lambda
    )

    # check dimensions are correct
    expect_equal(ncol(pr@matrix), 5)
    expect_equal(nrow(pr@matrix), 5)
  
})


test_that("pcor.bshrink is correct when ivar is alr",{

    message_test("pcor.bshrink is correct when ivar is alr")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    out <- propr:::basis_shrinkage(ct, outtype='alr')
    mat <- out$pcor
    lambda <- out$lambda
    attributes(mat) = NULL
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- suppressWarnings(propr(X, metric = "pcor.bshrink", ivar='alr'))

    # expect counts to have zeros replaced
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # NOTE that the data is not logratio transformed while computing pcor.bshrink
    # it is internally handled by covariance conversion.
    # so pr@logratio should be equal to pr@counts
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(pr@counts)
    )
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expect same lambda
    expect_equal(
      pr@lambda,
      lambda
    )

    # check dimensions are correct
    expect_equal(ncol(pr@matrix), 5)
    expect_equal(nrow(pr@matrix), 5)
  
})

test_that("test that pcor.bshrink gives error when ivar is NA", {

    message_test("test that pcor.bshrink gives error when ivar is NA")
    expect_error(
        propr(X, metric = "pcor.bshrink", ivar=NA)
    )
})

test_that("pcor.bshrink with alr and clr are the same", {

    message_test("pcor.bshrink with alr and clr are the same")

    # compute pcor with alr
    pr_alr <- suppressWarnings(propr(X, metric = "pcor.bshrink", ivar='alr'))
    pcor_alr <- getMatrix(pr_alr)[1:4, 1:4]

    # compute pcor with clr
    pr_clr <- propr(X, metric = "pcor.bshrink", ivar="clr")
    pcor_clr <- getMatrix(pr_clr)[1:4, 1:4]

    # expect that the coefficients are the same
    expect_equal(
        round(pcor_alr, 8),
        round(pcor_clr, 8)
    )

})
