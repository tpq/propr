library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# compute pcor.bshrink with alr
pr <- propr(X, metric = "pcor.bshrink", ivar="alr")
pcor_alr <- getMatrix(pr)[1:4, 1:4]

# compute pcor.bshrink with clr
pr <- propr(X, metric = "pcor.bshrink", ivar="clr")
pcor_clr <- getMatrix(pr)[1:4, 1:4]

# test that the results are the same
test_that("pcor.bshrink with alr and clr are the same", {

  expect_equal(
    pcor_alr,
    pcor_clr
  )
})
