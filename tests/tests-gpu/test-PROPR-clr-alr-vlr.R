library(testthat)
library(propr)

set.seed(101)

N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

matX <- as.matrix(X[])
tol <- 1e-6

test_that("propr:::clrRcpp CPU vs GPU", {
  cpu_clr <- tryCatch(
    propr:::clrRcpp(matX),
    error = function(e) skip("propr:::clrRcpp not available")
  )
  gpu_clr <- tryCatch(
    propr:::clrRcpp(matX, use_gpu = TRUE),
    error = function(e) skip("propr:::clrRcpp GPU variant not available")
  )
  expect_equal(cpu_clr, gpu_clr, tolerance = tol)
})

test_that("propr:::alrRcpp CPU vs GPU", {
  cpu_alr <- tryCatch(
    propr:::alrRcpp(matX, ivar = 5),
    error = function(e) skip("propr:::alrRcpp not available")
  )
  gpu_alr <- tryCatch(
    propr:::alrRcpp(matX, ivar = 5, use_gpu = TRUE),
    error = function(e) skip("propr:::alrRcpp GPU variant not available")
  )
  expect_equal(cpu_alr, gpu_alr, tolerance = tol)
})
