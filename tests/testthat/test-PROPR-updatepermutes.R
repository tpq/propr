library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# compute pcor.bshrink without fixed seed
pcorbshrink1 <- propr(X, metric = "pcor.bshrink", p=10)
pcorbshrink1 <- updateCutoffs(pcorbshrink1)

pcorbshrink2 <- propr(X, metric = "pcor.bshrink", p=10)
pcorbshrink2 <- updateCutoffs(pcorbshrink2)

set.seed(0)

# compute pcor.bshrink with fixed seed
pcorbshrink1_ <- propr(X, metric = "pcor.bshrink", p=10)
pcorbshrink1_ <- updateCutoffs(pcorbshrink1_)

set.seed(0)

pcorbshrink2_ <- propr(X, metric = "pcor.bshrink", p=10)
pcorbshrink2_ <- updateCutoffs(pcorbshrink2_)

# test that the results are as expected
test_that("test that fdr will stay the same when seed is given", {

  expect_false(
    isTRUE(all.equal(
        pcorbshrink1@permutes,
        pcorbshrink2@permutes
    ))
  )

  expect_false(
    isTRUE(all.equal(
        pcorbshrink1@fdr,
        pcorbshrink2@fdr
    ))
  )

  expect_equal(
    pcorbshrink1_@permutes,
    pcorbshrink2_@permutes
  )

  expect_equal(
    pcorbshrink1_@fdr,
    pcorbshrink2_@fdr
  )
})