library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

pr <- propr(X, metric = "rho")
rho <- getMatrix(pr)
diag(rho) <- 0

test_that("half-matrix correctly turned into matrix", {

  expect_equal(
    rho[1:16],
    propr:::half2mat(propr:::lltRcpp(rho))[1:16]
  )
})
