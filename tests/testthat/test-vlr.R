library(propr)
context("vlr")

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# VLR using X.clr
X.clr <- propr:::proprCLR(X)
Cov    <- stats::var(X.clr)
D      <- ncol(X.clr)
VarCol <- matrix(rep(diag(Cov), D), ncol = D)
VLR.clr <- -2 * Cov + VarCol + t(VarCol)

# VLR using X.alr
X.alr <- propr:::proprALR(X, ivar = 5)
Cov    <- stats::var(X.alr)
D      <- ncol(X.alr)
VarCol <- matrix(rep(diag(Cov), D), ncol = D)
VLR.alr <- -2 * Cov + VarCol + t(VarCol)

test_that("vlr shows subcompositional coherence", {

  expect_equivalent(
    propr:::proprVLR(X),
    propr:::proprVLR(X / rowSums(X))
  )

  expect_equivalent(
    propr:::proprVLR(X),
    VLR.clr
  )

  expect_equivalent(
    propr:::proprVLR(X)[-5, -5],
    VLR.alr
  )

  expect_equivalent(
    VLR.clr[-5, -5],
    VLR.alr
  )
})
