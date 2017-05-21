library(propr)
library(compositions)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

test_that("propr *lr matches compositions *lr", {

  expect_equivalent(
    round(
      unlist(propr:::proprCLR(X)),
      4),
    round(
      as.vector(compositions::clr(X)),
      4)
  )

  expect_equivalent(
    round(
      unlist(propr:::proprALR(X, ivar = 5)),
      4),
    round(
      as.vector(compositions::alr(X, ivar = 5)),
      4)
  )
})
