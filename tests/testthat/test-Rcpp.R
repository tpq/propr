library(propr)
context("Rcpp")

# Avoid modify-in-place behavior by using integer matrix
X <- t(data.frame("a" = sample(1:5), "b" = sample(1:5), "c" = sample(1:5),
                  "d" = sample(1:5), "e" = sample(1:5), "f" = sample(1:5)))

test_that("Rcpp backend matches original R backend", {

  expect_equivalent(
    cov(X),
    covRcpp(X)
  )

  expect_equivalent(
    propr:::proprVLR(X),
    vlrRcpp(X)
  )

  expect_equivalent(
    propr:::proprCLR(X),
    clrRcpp(X)
  )

  expect_error(
    alrRcpp(X, 0)
  )

  expect_equivalent(
    propr:::proprALR(X, 5),
    alrRcpp(X, 5)[, -5]
  )

  expect_equivalent(
    propr:::proprPhit(X),
    phiRcpp(X)
  )

  expect_equivalent(
    propr:::proprPerb(X, 0),
    rhoRcpp(X, clrRcpp(X[]), 0)
  )

  expect_equivalent(
    propr:::proprPerb(X, 5),
    rhoRcpp(X, alrRcpp(X[], 5), 5)[-5, -5]
  )
})

# Force modify-in-place behavior by using double matrix
X <- t(data.frame("a" = abs(rnorm(5)), "b" = abs(rnorm(5)), "c" = abs(rnorm(5)),
                  "d" = abs(rnorm(5)), "e" = abs(rnorm(5)), "f" = abs(rnorm(5))))

phi <- phit(X)
rho <- perb(X)

test_that("Rcpp backend does not distort phit or perb results", {

  expect_equivalent(
    phi@logratio,
    clrRcpp(X)
  )

  expect_equivalent(
    phi@logratio,
    X
  )

  expect_equivalent(
    rho@logratio,
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

mat <- phiRcpp(X)
tri <- propr:::proprTri(mat)

test_that("Rcpp backend indexes correct values", {

  expect_equivalent(
    tri,
    mat[indexPairs(mat, "all")]
  )

  expect_equivalent(
    tri[tri < 2],
    mat[indexPairs(mat, "<", 2)]
  )

  expect_equivalent(
    tri[tri > 2],
    mat[indexPairs(mat, ">", 2)]
  )
})
