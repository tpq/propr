library(propr)

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

# Avoid modify-in-place behavior by using integer matrix
X <- t(data.frame("a" = sample(1:5), "b" = sample(1:5), "c" = sample(1:5),
                  "d" = sample(1:5), "e" = sample(1:5), "f" = sample(1:5)))

test_that("Rcpp backend matches original R backend", {

  expect_equivalent(
    cov(X),
    propr:::covRcpp(X)
  )

  expect_equivalent(
    propr:::proprVLR(X),
    propr:::vlrRcpp(X)
  )

  expect_equivalent(
    propr:::proprCLR(X),
    propr:::clrRcpp(X)
  )

  expect_error(
    propr:::alrRcpp(X, 0)
  )

  expect_equivalent(
    propr:::proprALR(X, 5),
    propr:::alrRcpp(X, 5)[, -5]
  )

  expect_equivalent(
    propr:::proprPhit(X),
    propr:::phiRcpp(X)
  )

  expect_equivalent(
    propr:::proprPerb(X, 0),
    propr:::rhoRcpp(X, propr:::clrRcpp(X[]), 0)
  )

  expect_equivalent(
    propr:::proprPerb(X, 5),
    propr:::rhoRcpp(X, propr:::alrRcpp(X[], 5), 5)[-5, -5]
  )
})

# Force modify-in-place behavior by using double matrix
X <- t(data.frame("a" = abs(rnorm(5)), "b" = abs(rnorm(5)), "c" = abs(rnorm(5)),
                  "d" = abs(rnorm(5)), "e" = abs(rnorm(5)), "f" = abs(rnorm(5))))

phi <- phit(X)
rho <- perb(X)

test_that("Rcpp backend does not distort phit or perb results", {

  expect_equivalent(
    as.matrix(phi@logratio),
    propr:::clrRcpp(X)
  )

  expect_equivalent(
    as.matrix(phi@logratio),
    X
  )

  expect_equivalent(
    as.matrix(rho@logratio),
    X
  )

  expect_equivalent(
    phi@counts,
    rho@counts
  )
})

# Force modify-in-place behavior by using double matrix
X <- t(data.frame("a" = abs(rnorm(5)), "b" = abs(rnorm(5)), "c" = abs(rnorm(5)),
                  "d" = abs(rnorm(5)), "e" = abs(rnorm(5)), "f" = abs(rnorm(5))))

mat <- propr:::phiRcpp(X)
tri <- propr:::proprTri(mat)

test_that("Rcpp backend indexes correct values", {

  expect_equivalent(
    tri,
    mat[propr:::indexPairs(mat, "all")]
  )

  expect_equivalent(
    tri[tri < 2],
    mat[propr:::indexPairs(mat, "<", 2)]
  )

  expect_equivalent(
    tri[tri > 2],
    mat[propr:::indexPairs(mat, ">", 2)]
  )
})
