library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# get phi
pr <- propr(X, metric = "phi")
phi <- getMatrix(pr)

# get rho
pr <- propr(X, metric = "rho")
rho <- getMatrix(pr)

# get phs
pr <- propr(X, metric = "phs")
phs <- getMatrix(pr)

# get cor
pr <- propr(X, metric = "cor")
cor <- getMatrix(pr)

# get beta
counts.clr <- t(apply(X, 1, function(x) log(x) - mean(log(x))))
counts.clr.var <- apply(counts.clr, 2, var)
A_j <- matrix(rep(counts.clr.var, length(counts.clr.var)), nrow = length(counts.clr.var))
A_i <- counts.clr.var
beta <- sqrt(sweep(A_j, 2, A_i, "/"))

# calculate alt rho, phs, cor
rho_ <- 1 - phi / (1 + beta^2)
phs_ <- (1 - rho_) / (1 + rho_)
cor_ <- (1 + beta^2 - phi) / (2 * beta)

test_that("calculating rho from phi matches rho", {

  expect_equal(
    rho_,
    rho
  )

  expect_equal(
    phs_,
    phs
  )

  expect_equal(
    cor_,
    cor
  )
})
