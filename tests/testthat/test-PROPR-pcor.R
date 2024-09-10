library(testthat)
library(propr)

# define data
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)
for (i in 1:4){
  X[sample(1:N, 10),i] <- 0
}

test_that("pcor is correct when ivar is given", {

  d1 <- list("clr", 5, c(1,3))
  d2 <- list(c(1:5), 5, c(1,3))

  for (i in 1:length(d1)){

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    lr <- logratio_without_alpha(ct, d2[[i]])
    cov <- cov(lr)
    mat <- suppressWarnings(corpcor::cor2pcor(cov))

    # compute pcor
    pr <- suppressWarnings(propr(X, metric = "pcor", ivar=d1[[i]]))

    # expect counts to have zeros replaced
    expect_true(
      all(pr@counts == ct)
    )

    # expect that logratio are equal
    expect_true(
      all(pr@logratio == lr)
    )

    # expect computed coefficients are equal
    expect_true(
      all(round(pr@matrix, 8) == round(mat, 8), na.rm=T)
    )

    # check shrinkage is not applied
    expect_equal(
      pr@lambda,
      NULL
    )
  }
})


test_that("pcor is correct when ivar is NA using previously transformed data", {

  d1 <- list("clr", 5, c(1,3))
  d2 <- list(c(1:5), 5, c(1,3))

  for (i in 1:length(d1)){

    # compute pcor manually
    ct <- simple_zero_replacement(X)
    lr <- logratio_without_alpha(ct, d2[[i]])
    cov <- cov(lr)
    mat <- suppressWarnings(corpcor::cor2pcor(cov))

    # compute pcor
    pr <- suppressWarnings(propr(lr, metric = "pcor", ivar=NA))

    # expect counts to contain the previously transformed data
    expect_true(
      all(pr@counts == lr)
    )

    # expect that the logratio slot also contains the exact input data
    expect_true(
      all(pr@logratio == lr)
    )

    # expect computed coefficients are equal
    expect_true(
      all(round(pr@matrix, 8) == round(mat, 8), na.rm=T)
    )

    # check shrinkage is not applied
    expect_equal(
      pr@lambda,
      NULL
    )
  }
})

test_that("pcor with alr and clr are the same", {

    # compute pcor with alr
    pr <- suppressWarnings(propr(X, metric = "pcor", ivar=5))
    pcor_alr <- getMatrix(pr)[1:4, 1:4]

    # compute pcor with clr
    pr <- propr(X, metric = "pcor", ivar="clr")
    pcor_clr <- getMatrix(pr)[1:4, 1:4]

    expect_equal(
      round(pcor_alr, 8),
      round(pcor_clr, 8)
    )

})
