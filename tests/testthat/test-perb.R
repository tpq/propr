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

# Calculate beta
counts.clr <- propr:::proprCLR(X)
counts.clr.var <- apply(counts.clr, 2, var)
A_j <- matrix(rep(counts.clr.var, length(counts.clr.var)), nrow = length(counts.clr.var))
A_i <- counts.clr.var
beta <- sqrt(sweep(A_j, 2, A_i, "/"))

# Calculate rho
rho <- 1 - phi/(1 + beta^2)

test_that("calculating rho from phi matches propr:::proprPerb", {

  expect_equal(
    rho,
    propr:::proprPerb(X, ivar = 0)
  )
})

test_that("calculating rho from perb matches propr:::proprPerb", {

  expect_equal(
    as.vector(perb(X)@matrix),
    as.vector(propr:::proprPerb(X, ivar = 0))
  )

  expect_equal(
    as.vector(perb(X, 3)@matrix[c(1:2, 4:5), c(1:2, 4:5)]),
    as.vector(propr:::proprPerb(X, ivar = 3))
  )
})

test_that("perb accepts ivar name or index", {

  expect_equal(
    perb(X, 2)@matrix,
    perb(X, "b")@matrix
  )
})

test_that("perb select argument works (clr)", {

  include <- c("b", "d", "e")
  mat1 <- perb(X, select = include)
  mat2 <- perb(X)
  mat2 <- subset(mat2, select = include)

  expect_equal(
    mat1@counts,
    mat2@counts
  )

  expect_equal(
    mat1@logratio,
    mat2@logratio
  )

  expect_equal(
    mat1@matrix,
    mat2@matrix
  )
})

test_that("perb select argument works (alr) without dropping ivar", {

  include <- c("b", "d", "e")
  mat1 <- perb(X, 4, select = include)
  mat2 <- perb(X, 4)
  mat2 <- subset(mat2, select = include)

  expect_equal(
    mat1@counts,
    mat2@counts
  )

  expect_equal(
    mat1@logratio,
    mat2@logratio
  )

  expect_equal(
    mat1@matrix,
    mat2@matrix
  )
})

test_that("perb select argument works (alr) with dropping ivar", {

  include <- c("b", "d", "e")
  mat1 <- perb(X, 3, select = include)
  mat2 <- perb(X, 3)
  mat2 <- subset(mat2, select = include)

  expect_equal(
    mat1@counts,
    mat2@counts
  )

  expect_equal(
    mat1@logratio,
    mat2@logratio
  )

  expect_equal(
    mat1@matrix,
    mat2@matrix
  )
})

test_that("perb select argument works (alr) out of order", {

  include <- c("e", "d", "b")
  mat1 <- perb(X, 3, select = include)
  mat2 <- perb(X, 3)
  mat2 <- subset(mat2, select = include)

  expect_equal(
    mat1@counts,
    mat2@counts
  )

  expect_equal(
    mat1@logratio,
    mat2@logratio
  )

  expect_equal(
    mat1@matrix,
    mat2@matrix
  )
})
