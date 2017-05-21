library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# Calculate phi
phi <- propr:::proprPhit(X, symmetrize = FALSE)

test_that("calculating phi from phit matches propr:::proprPhit", {

  expect_equal(
    as.vector(phit(X, symmetrize = FALSE)@matrix),
    as.vector(phi)
  )
})

# Calculate r
counts.clr <- propr:::proprCLR(X)
counts.clr.var <- apply(counts.clr, 2, var)
A_j <- matrix(rep(counts.clr.var, length(counts.clr.var)), nrow = length(counts.clr.var))
A_i <- counts.clr.var
beta <- sqrt(sweep(A_j, 2, A_i, "/"))
r <- (1 + beta^2 - phi) / (2 * beta)

test_that("calculating r from phi matches stats::cor", {

  expect_equal(
    r,
    stats::cor(propr:::proprCLR(X))
  )
})
