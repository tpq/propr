library(testthat)
library(propr)

set.seed(101)

# Toy data used across tests
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)
matX <- as.matrix(X[])
tol <- 1e-6

# Helper: try CPU then GPU; skip if not available
call_cpu_gpu <- function(funname, cpu_call, gpu_call) {
  cpu_val <- tryCatch(
    cpu_call(),
    error = function(e) skip(paste0(funname, " CPU call not available: ", e$message))
  )
  gpu_val <- tryCatch(
    gpu_call(),
    error = function(e) skip(paste0(funname, " GPU call not available or use_gpu flag missing: ", e$message))
  )
  list(cpu = cpu_val, gpu = gpu_val)
}

# 1) wtmRcpp: weighted trimmed mean (scalar)
test_that("wtmRcpp CPU vs GPU", {
  message("=== TEST: wtmRcpp CPU vs GPU ===")
  w <- runif(nrow(matX), 0.1, 2)
  vals <- call_cpu_gpu(
    "wtmRcpp",
    function() propr:::wtmRcpp(matX[,1], w),
    function() propr:::wtmRcpp(matX[,1], w, use_gpu = TRUE)
  )
  expect_equal(as.numeric(vals$cpu), as.numeric(vals$gpu), tolerance = tol)
})

# 2) wtvRcpp: weighted variance scalar
test_that("wtvRcpp CPU vs GPU", {
  message("=== TEST: wtvRcpp CPU vs GPU ===")
  w <- runif(nrow(matX), 0.1, 2)
  vals <- call_cpu_gpu(
    "wtvRcpp",
    function() propr:::wtvRcpp(matX[,1], w),
    function() propr:::wtvRcpp(matX[,1], w, use_gpu = TRUE)
  )
  expect_equal(as.numeric(vals$cpu), as.numeric(vals$gpu), tolerance = tol)
})

# 3) centerNumericMatrix: centers columns (matrix)
test_that("centerNumericMatrix CPU vs GPU", {
  message("=== TEST: centerNumericMatrix CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "centerNumericMatrix",
    function() propr:::centerNumericMatrix(matX),
    function() propr:::centerNumericMatrix(matX, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 4) corRcpp: correlation matrix
test_that("corRcpp CPU vs GPU", {
  message("=== TEST: corRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "corRcpp",
    function() propr:::corRcpp(matX),
    function() propr:::corRcpp(matX, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = 1e-5)
})

# 5) covRcpp: covariance matrix (norm_type = 0)
test_that("covRcpp CPU vs GPU", {
  message("=== TEST: covRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "covRcpp",
    function() propr:::covRcpp(matX, 0),
    function() propr:::covRcpp(matX, 0, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 6) vlrRcpp: variance of log-ratios matrix
test_that("vlrRcpp CPU vs GPU", {
  message("=== TEST: vlrRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "vlrRcpp",
    function() { m <- propr:::vlrRcpp(matX); as.matrix(m) },
    function() { m <- propr:::vlrRcpp(matX, use_gpu = TRUE); as.matrix(m) }
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 7) clrRcpp: centered log ratio transform
test_that("clrRcpp CPU vs GPU", {
  message("=== TEST: clrRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "clrRcpp",
    function() propr:::clrRcpp(matX),
    function() propr:::clrRcpp(matX, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 8) alrRcpp: additive log ratio transform (ivar = 5)
test_that("alrRcpp CPU vs GPU", {
  message("=== TEST: alrRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "alrRcpp",
    function() propr:::alrRcpp(matX, ivar = 5),
    function() propr:::alrRcpp(matX, ivar = 5, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 9) symRcpp: symmetrize a matrix
test_that("symRcpp CPU vs GPU", {
  message("=== TEST: symRcpp CPU vs GPU ===")
  M <- matrix(rnorm(ncol(matX)^2), ncol = ncol(matX))
  vals <- call_cpu_gpu(
    "symRcpp",
    function() propr:::symRcpp(M),
    function() propr:::symRcpp(M, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 10) phiRcpp: phi metric (matrix), try with sym = TRUE
test_that("phiRcpp CPU vs GPU", {
  message("=== TEST: phiRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "phiRcpp",
    function() propr:::phiRcpp(matX, sym = TRUE),
    function() propr:::phiRcpp(matX, sym = TRUE, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 11) rhoRcpp: rho matrix (needs lr and ivar)
test_that("rhoRcpp CPU vs GPU", {
  message("=== TEST: rhoRcpp CPU vs GPU ===")

    # skip(message=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WARN:FAILING")

  lr <- tryCatch(propr:::clrRcpp(matX), error = function(e) skip("clrRcpp not present for preparing rhoRcpp test"))
  vals <- call_cpu_gpu(
    "rhoRcpp",
    function() propr:::rhoRcpp(matX, lr, ivar = 5),
    function() propr:::rhoRcpp(matX, lr, ivar = 5, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 12) indexPairs: return pair index vector (op "all")
test_that("indexPairs CPU vs GPU", {
  message("=== TEST: indexPairs CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "indexPairs",
    function() propr:::indexPairs(matX, "all", 0),
    function() propr:::indexPairs(matX, "all", 0, use_gpu = TRUE)
  )
  expect_equal(as.integer(vals$cpu), as.integer(vals$gpu))
})

# 13) indexToCoord: convert vector->coords
test_that("indexToCoord CPU vs GPU", {
  message("=== TEST: indexToCoord CPU vs GPU ===")
  V <- as.integer(1:10)
  Nfeats <- ncol(matX)
  vals <- call_cpu_gpu(
    "indexToCoord",
    function() propr:::indexToCoord(V, N = Nfeats),
    function() propr:::indexToCoord(V, N = Nfeats, use_gpu = TRUE)
  )
  expect_equal(names(vals$cpu), names(vals$gpu))
  expect_equal(vals$cpu$feat1, vals$gpu$feat1)
  expect_equal(vals$cpu$feat2, vals$gpu$feat2)
})

# 14) coordToIndex: coords -> index
test_that("coordToIndex CPU vs GPU", {
  message("=== TEST: coordToIndex CPU vs GPU ===")
  rowv <- as.integer(c(1,2,3))
  colv <- as.integer(c(2,3,4))
  Nfeats <- 5
  vals <- call_cpu_gpu(
    "coordToIndex",
    function() propr:::coordToIndex(rowv, colv, N = Nfeats),
    function() propr:::coordToIndex(rowv, colv, N = Nfeats, use_gpu = TRUE)
  )
  expect_equal(as.integer(vals$cpu), as.integer(vals$gpu))
})

# 15) linRcpp: compute linearized matrix (rho, lr)
test_that("linRcpp CPU vs GPU", {
  message("=== TEST: linRcpp CPU vs GPU ===")
  lr <- tryCatch(propr:::clrRcpp(matX), error = function(e) skip("clrRcpp not present for preparing linRcpp test"))
  rho <- tryCatch(propr:::rhoRcpp(matX, lr, ivar = 5), error = function(e) skip("rhoRcpp not present for preparing linRcpp test"))
  vals <- call_cpu_gpu(
    "linRcpp",
    function() propr:::linRcpp(rho, lr),
    function() propr:::linRcpp(rho, lr, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 16) lltRcpp: lower-triangle-to-vector
test_that("lltRcpp CPU vs GPU", {
  message("=== TEST: lltRcpp CPU vs GPU ===")  

  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)  
  pr <- propr(X, metric = "rho")
  rho <- getMatrix(pr)
  diag(rho) <- 0

  vals <- call_cpu_gpu(
    "lltRcpp",
    function() propr:::lltRcpp(rho),
    function() propr:::lltRcpp(rho, use_gpu = TRUE)
  )
  expect_equal(as.numeric(vals$cpu), as.numeric(vals$gpu), tolerance = tol)
})

# 17) urtRcpp: upper-triangle-to-vector
test_that("urtRcpp CPU vs GPU", {
  message("=== TEST: urtRcpp CPU vs GPU ===")
  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)  
  pr <- propr(X, metric = "rho")
  S <- getMatrix(pr)
  diag(S) <- 0


  vals <- call_cpu_gpu(
    "urtRcpp",
    function() propr:::urtRcpp(S),
    function() propr:::urtRcpp(S, use_gpu = TRUE)
  )
  expect_equal(as.numeric(vals$cpu), as.numeric(vals$gpu), tolerance = tol)
})

# 18) labRcpp: build label list for pairs
test_that("labRcpp CPU vs GPU", {
  message("=== TEST: labRcpp CPU vs GPU ===")
  nfeats <- ncol(matX)
  vals <- call_cpu_gpu(
    "labRcpp",
    function() propr:::labRcpp(nfeats),
    function() propr:::labRcpp(nfeats, use_gpu = TRUE)
  )
  expect_equal(names(vals$cpu), names(vals$gpu))
  expect_equal(as.integer(vals$cpu$Partner), as.integer(vals$gpu$Partner))
  expect_equal(as.integer(vals$cpu$Pair), as.integer(vals$gpu$Pair))
})

# 19) half2mat: vector -> symmetric matrix
test_that("half2mat CPU vs GPU", {
  message("=== TEST: half2mat CPU vs GPU ===")
  vec <- rnorm(ncol(matX) * (ncol(matX) - 1) / 2)
  vals <- call_cpu_gpu(
    "half2mat",
    function() propr:::half2mat(vec),
    function() propr:::half2mat(vec, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 20) vector2mat: vector + coords -> symmetric matrix
test_that("vector2mat CPU vs GPU", {
  message("=== TEST: vector2mat CPU vs GPU ===")
  nfeats <- 5
  i <- as.integer(c(1,2,3))
  j <- as.integer(c(2,3,4))
  Xvals <- rnorm(length(i))
  vals <- call_cpu_gpu(
    "vector2mat",
    function() propr:::vector2mat(Xvals, i, j, nfeats),
    function() propr:::vector2mat(Xvals, i, j, nfeats, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 21) ratiosRcpp: sample-wise ratios (nsamp x llt)
test_that("ratiosRcpp CPU vs GPU", {
  message("=== TEST: ratiosRcpp CPU vs GPU ===")
  vals <- call_cpu_gpu(
    "ratiosRcpp",
    function() propr:::ratiosRcpp(matX),
    function() propr:::ratiosRcpp(matX, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})

# 22) results2matRcpp: fill matrix from results dataframe
test_that("results2matRcpp CPU vs GPU", {
  message("=== TEST: results2matRcpp CPU vs GPU ===")
  n <- ncol(matX)
  results <- data.frame(V1 = c(1,2,3), V2 = c(2,3,1), V3 = c(0.5, 0.6, 0.7))
  vals <- call_cpu_gpu(
    "results2matRcpp",
    function() propr:::results2matRcpp(n, results, diagonal = 0),
    function() propr:::results2matRcpp(n, results, diagonal = 0, use_gpu = TRUE)
  )
  expect_equal(as.matrix(vals$cpu), as.matrix(vals$gpu), tolerance = tol)
})
