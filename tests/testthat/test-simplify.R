library(propr)

set.seed(1)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * abs(rnorm(N, mean = 1, sd = 10))
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# Index and simplify
phi <- phit(X)
phi.index <- phi["<", .5]
simple <- simplify(phi.index)

test_that("simplify correctly re-indexes proportionality matrix", {

  expect_equal(
    simple@matrix[simple@pairs],
    phi.index@matrix[phi.index@pairs]
  )
})
