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

test_that("pcor.shrink is correct when ivar is clr",{

    message_test("pcor.shrink is correct when ivar is clr")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    clr <- t(apply(ct, 1, function(x) log(x) - mean(log(x))))
    mat <- corpcor::pcor.shrink(clr)
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- propr(X, metric = "pcor.shrink", ivar="clr")

    # expect counts to have zeros replaced
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # expect that clr are equal
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(clr)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )
  
})

test_that("pcor.shrink gives error when ivar is alr", {
  
    message_test("pcor.shrink gives error when ivar is alr")
    expect_error(
      propr(X, metric = "pcor.shrink", ivar='alr')
    )
})

test_that("pcor.shrink is correct when ivar is 5",{

    message_test("pcor.shrink is correct when ivar is 5")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    alr <- t(apply(ct, 1, function(x) log(x) - log(x[5])))
    mat <- suppressWarnings(corpcor::pcor.shrink(alr))
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- suppressWarnings(propr(X, metric = "pcor.shrink", ivar=5))

    # expect counts to have zeros replaced
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # expect that alr are equal
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(alr)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )
  
})

test_that("pcor.shrink is correct when ivar is 1,3",{

    message_test("pcor.shrink is correct when ivar is 1,3")
  
    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    alr <- logratio_without_alpha(ct, c(1,3))
    mat <- suppressWarnings(corpcor::pcor.shrink(alr))
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- suppressWarnings(propr(X, metric = "pcor.shrink", ivar=c(1,3)))

    # expect counts to have zeros replaced
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # expect that alr are equal
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(alr)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )
  
})

test_that("pcor.shrink is correct when ivar is NA using previously clr transformed data", {

    message_test("pcor.shrink is correct when ivar is NA using previously clr transformed data")

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    clr <- t(apply(ct, 1, function(x) log(x) - mean(log(x))))
    mat <- corpcor::pcor.shrink(clr)
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- propr(clr, metric = "pcor.shrink", ivar=NA)

    # expect counts to contain the previously transformed data
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(clr)
    )

    # expect that the logratio slot also contains the exact input data
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(clr)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )

})

test_that("pcor.shrink is correct when ivar is NA using previously alr transformed data", {

    message_test("pcor.shrink is correct when ivar is NA using previously alr transformed data")

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    alr <- t(apply(ct, 1, function(x) log(x) - log(x[5])))
    mat <- suppressWarnings(corpcor::pcor.shrink(alr))
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- suppressWarnings(propr(alr, metric = "pcor.shrink", ivar=NA))

    # expect counts to contain the previously transformed data
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(alr)
    )

    # expect that the logratio slot also contains the exact input data
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(alr)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )

})

test_that("pcor.shrink is correct when ivar is NA using raw data", {

    message_test("pcor.shrink is correct when ivar is NA using raw data")

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    rownames(ct) <- rownames(X)
    mat <- suppressWarnings(corpcor::pcor.shrink(ct))
    lambda <- attr(mat, "lambda")
    attributes(mat) = NULL
    mat <- matrix(mat, ncol=ncol(X), nrow=ncol(X))
    class(mat) <- "matrix"
    colnames(mat)  <- colnames(X)
    rownames(mat)  <- colnames(X)

    # compute pcor
    pr <- suppressWarnings(propr(ct, metric = "pcor.shrink", ivar=NA))

    # expect counts to contain the raw data
    expect_equal(
      as.matrix(pr@counts),
      as.matrix(ct)
    )

    # expect that the logratio slot also contains the exact input data
    expect_equal(
      as.matrix(pr@logratio),
      as.matrix(ct)
    )   
  
    # expect computed coefficients are equal
    expect_equal(
      round(pr@matrix, 8),
      round(mat, 8)
    )

    # expected same shrinkage lambda
    expect_equal(
      pr@lambda,
      lambda
    )

})

test_that("pcor.shrink with alr and clr are not the same because of different shrinkage estimation", {

    message_test("pcor.shrink with alr and clr are not the same because of different shrinkage estimation")

    # compute pcor with alr
    pr_alr <- suppressWarnings(propr(X, metric = "pcor.shrink", ivar=5))
    pcor_alr <- getMatrix(pr_alr)[1:4, 1:4]

    # compute pcor with clr
    pr_clr <- propr(X, metric = "pcor.shrink", ivar="clr")
    pcor_clr <- getMatrix(pr_clr)[1:4, 1:4]

    # expect that the coefficients are not equal
    expect_false(identical(
        round(pcor_alr, 8),
        round(pcor_clr, 8)
    ))

    # expect that the shrinkage lambda are not equal
    expect_false(identical(
        pr_alr@lambda,
        pr_clr@lambda
    ))
})
