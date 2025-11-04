library(testthat)
library(propr)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)
tol <- 1e-6

# VLR using propr with CLR ivar
X.clr <- propr(X, ivar = "clr")@logratio
Cov    <- stats::var(X.clr)
D      <- ncol(X.clr)
VarCol <- matrix(rep(diag(Cov), D), ncol = D)
VLR.clr <- -2 * Cov + VarCol + t(VarCol)

# VLR using propr with ALR ivar
X.alr <- propr(X, ivar = 5)@logratio
Cov    <- stats::var(X.alr)
D      <- ncol(X.alr)
VarCol <- matrix(rep(diag(Cov), D), ncol = D)
VLR.alr <- -2 * Cov + VarCol + t(VarCol)

# VLR using propr VLR functions
X.vlr1 <- propr(X, metric = "vlr")@matrix
X.vlr2 <- propr(X / rowSums(X), metric = "vlr")@matrix
X.vlr3 <- propr:::vlrRcpp(as.matrix(X[]))
rownames(X.vlr3) <- colnames(X)
colnames(X.vlr3) <- colnames(X)
X.vlr4 <- propr:::lr2vlr(as.matrix(X.clr))
rownames(X.vlr4) <- colnames(X)
colnames(X.vlr4) <- colnames(X)

test_that("vlr shows subcompositional coherence", {

  expect_equal(
    X.vlr1, # propr(how = "vlr")
    X.vlr2 # propr(how = "vlr")
  )

  expect_equal(
    X.vlr1, # propr(how = "vlr")
    X.vlr3 # vlrRcpp
  )

  expect_equal(
    X.vlr3, # vlrRcpp
    X.vlr4 # lr2vlr
  )

  expect_equal(
    X.vlr3, # vlrRcpp
    VLR.clr
  )

  expect_equal(
    X.vlr3, # vlrRcpp
    VLR.alr
  )

  expect_equal(
    VLR.clr,
    VLR.alr
  )
})

test_that("propr:::clrRcpp CPU vs GPU", {
  cpu_clr <- tryCatch(
    propr:::clrRcpp(as.matrix(X[])),
    error = function(e) skip(paste0("propr:::clrRcpp not available — ", conditionMessage(e)))
  )
  gpu_clr <- tryCatch(
    propr:::clrRcpp(as.matrix(X[]), use_gpu = TRUE),
    error = function(e) skip(paste0("propr:::clrRcpp GPU variant not available — ", conditionMessage(e)))
  )
  expect_equal(cpu_clr, gpu_clr, tolerance = tol)
})

test_that("propr:::alrRcpp CPU vs GPU", {
  cpu_alr <- tryCatch(
    propr:::alrRcpp(as.matrix(X[]), ivar = 5),
    error = function(e) skip(paste0("propr:::alrRcpp not available — ", conditionMessage(e)))
  )
  gpu_alr <- tryCatch(
    propr:::alrRcpp(as.matrix(X[]), ivar = 5, use_gpu = TRUE),
    error = function(e) skip(paste0("propr:::alrRcpp GPU variant not available — ", conditionMessage(e)))
  )
  expect_equal(cpu_alr, gpu_alr, tolerance = tol)
})
