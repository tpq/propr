library(testthat)
library(propr)
library(Rcpp)

# Create matrices of 0 and 1
A <- matrix(c(1, 1, 0, 1, 0, 
              1, 1, 1, 0, 1, 
              0, 1, 1, 0, 1, 
              1, 0, 0, 1, 1, 
              0, 1, 1, 1, 1), 
            nrow = 5, byrow = TRUE)
K <- matrix(c(1, 0, 
              0, 1, 
              1, 0, 
              0, 1, 
              1, 1), 
            nrow = 5, byrow = TRUE)
colnames(K) <- c("C1", "C2")
rownames(A) <- c('a', 'b', 'c', 'd', 'e')
rownames(K) <- c('a', 'b', 'c', 'd', 'e')

# ===================== #
# old code for graflex  #
# ===================== #

cppFunction("
  IntegerVector sample_idx(int n) {
    IntegerVector idx = sample(n, n, false) - 1;
    return idx;
  }
")

permuteOR_old <- function(A, G, p = 500) {
  Gstar <- G[lower.tri(G)]
  res <- lapply(1:p, function(i) {
    # index <- sample(ncol(A), ncol(A), replace = FALSE)
    index <- sample_idx(ncol(A))+1  # use the same random sampling generator
    A <- A[index, index]
    Astar <- A[lower.tri(A)]
    getOR_old(Astar, Gstar)
  })

  do.call("rbind", res)
}

binTab_old <- function(A, G) {
  diff <- A != G
  only1 <- A[diff]
  b <- as.numeric(sum(only1))
  c <- as.numeric(length(only1) - b)

  same <- !diff
  double1 <- A[same]
  a <- as.numeric(sum(double1))
  d <- as.numeric(length(double1) - a)

  matrix(c(d, b, c, a), 2, 2)
}

getOR_old <- function(A, G) {
  tab <- binTab_old(A, G)
  or <- (tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1])
  data.frame(
    "Neither" = tab[1, 1],
    "G.only" = tab[1, 2],
    "A.only" = tab[2, 1],
    "Both" = tab[2, 2],
    "Odds" = or,
    "LogOR" = log(or)
  )
}

calculateOR_old <- function(A, G) {
  Astar <- A[lower.tri(A)]
  Gstar <- G[lower.tri(G)]
  getOR_old(Astar, Gstar)
}

getFDR_old <- function(actual, permuted) {
  actual$FDR.under <-
    sum(permuted$Odds <= actual$Odds) / nrow(permuted)
  actual$FDR.over <-
    sum(permuted$Odds >= actual$Odds) / nrow(permuted)
  actual
}

runGraflex_old <- function(A, K, p = 500) {
  if (nrow(A) != nrow(K))
    stop("'A' and 'K' must have identical rows.")

  numTicks <- 0
  res <- lapply(1:ncol(K), function(k) {
    Gk <- K[, k] %*% t(K[, k])
    actual <- calculateOR_old(A, Gk)
    permuted <- permuteOR_old(A, Gk, p = p)
    actual <- getFDR_old(actual, permuted)
    actual$Permutes <- p
    actual$Concept <- colnames(K)[k]
    actual
  })

  do.call("rbind", res)
}

# ===================== #
# run tests             #
# ===================== #

test_that("check if odds ratio is computed correctly", {

  G <- K[, 1] %*% t(K[, 1])
  A1 <- A[lower.tri(A)]
  G1 <- G[lower.tri(G)]

  # Expected values
  a <- length(A1[A1 == 0 & G1 == 0])  # not in A not in G
  b <- length(A1[A1 == 0 & G1 == 1])  # not in A but in G
  c <- length(A1[A1 == 1 & G1 == 0])  # in A but not in G
  d <- length(A1[A1 == 1 & G1 == 1])  # in A in G
  expected_odds_ratio <- (a * d) / (b * c)

  # compute odds ratio table
  or_table <- propr:::getOR(A, G)

  # check
  expect_equal(a, as.integer(or_table[1]))
  expect_equal(b, as.integer(or_table[2]))
  expect_equal(c, as.integer(or_table[3]))
  expect_equal(d, as.integer(or_table[4]))
  expect_equal(expected_odds_ratio, or_table[5])
})

test_that("check that getFDR works properly", {

  # create a vector of values between 0 and 1
  x <- runif(100)
  actual <- runif(1)

  # compute FDR
  res <- propr:::getFDR(actual, x)

  # check
  expect_equal(res$under, sum(x <= actual) / length(x))
  expect_equal(res$over, sum(x >= actual) / length(x))
})


test_that("check that getFDR works properly GPU", {

  # create a vector of values between 0 and 1
  x <- runif(100)
  actual <- runif(1)

  # compute FDR
  res <- propr:::getFDR(actual, x, use_gpu=TRUE)

  # check
  expect_equal(res$under, sum(x <= actual) / length(x))
  expect_equal(res$over, sum(x >= actual) / length(x))
})

test_that("check if runGraflex produce the expected contingency table and odds ratio values", {
  
  Astar <- A[lower.tri(A)]

  # Expected values for C1
  G1 <- K[, 1] %*% t(K[, 1])
  G1 <- G1[lower.tri(G1)]
  a1 <- length(which(Astar == 0 & G1 == 0))  # not in A not in G
  b1 <- length(which(Astar == 0 & G1 == 1))  # not in A but in G
  c1 <- length(which(Astar == 1 & G1 == 0))  # in A but not in G
  d1 <- length(which(Astar == 1 & G1 == 1))  # in A in G
  expected_odds_ratio1 <- (a1 * d1) / (b1 * c1)

  # Expected values for C2
  G2 <- K[, 2] %*% t(K[, 2])
  G2 <- G2[lower.tri(G2)]
  a2 <- length(which(Astar == 0 & G2 == 0))  # not in A not in G
  b2 <- length(which(Astar == 0 & G2 == 1))  # not in A but in G
  c2 <- length(which(Astar == 1 & G2 == 0))  # in A but not in G
  d2 <- length(which(Astar == 1 & G2 == 1))  # in A in G
  expected_odds_ratio2 <- (a2 * d2) / (b2 * c2)

  # compute graflex
  res <- propr:::runGraflex(A, K)

  # check
  expect_equal(a1, as.numeric(res[1,1]))
  expect_equal(b1, as.numeric(res[1,2]))
  expect_equal(c1, as.numeric(res[1,3]))
  expect_equal(d1, as.numeric(res[1,4]))
  expect_equal(expected_odds_ratio1, as.numeric(res[1,5]))
  expect_equal(a2, as.numeric(res[2,1]))
  expect_equal(b2, as.numeric(res[2,2]))
  expect_equal(c2, as.numeric(res[2,3]))
  expect_equal(d2, as.numeric(res[2,4]))
  expect_equal(expected_odds_ratio2, as.numeric(res[2,5]))
})

test_that("check reproducibility seed works", {

    set.seed(0)
    res1 <- propr:::runGraflex(A, K, p=100)
    set.seed(0)
    res2 <- propr:::runGraflex(A, K, p=100)
    set.seed(123)
    res3 <- propr:::runGraflex(A, K, p=100)
  
    expect_equal(res1, res2)
    expect_false(isTRUE(all.equal(res1, res3)))
})

test_that("check that permuteOR works as the old code", {

  G <- K[, 1] %*% t(K[, 1])
  
  # compute permuted odds ratios
  set.seed(0)
  res1 <- propr:::permuteOR(A, G, 10, use_gpu=TRUE)
  set.seed(0)
  res2 <- permuteOR_old(A, G, 10)

  df1 <- as.data.frame(res1[, 1:6])
  colnames(df1) <- c("Neither", "G.only", "A.only", "Both", "Odds", "LogOR")
    
  # check
  expect_true(all(res1[,-c(7,8)] == res2))
})

test_that("check if same results are obtained compared to old code", {

    # compute graflex
    set.seed(0)
    res <- propr:::runGraflex(A, K, p=100)
    res$LogOR <- round(res$LogOR, 5)
    
    # compute graflex with old code
    set.seed(0)
    res_old <- runGraflex_old(A, K, p=100)
    res_old$LogOR <- round(res_old$LogOR, 5)
    
    # check
    expect_equal(res, res_old)
})

test_that("check that when ncores > 1 works", {

  pp <- 100

  # compute graflex
  res1 <- propr:::runGraflex(A, K, p=pp, ncores=1)
  res2 <- propr:::runGraflex(A, K, p=pp, ncores=2)

  # check
  cols <- c("Neither", "G.only", "A.only", "Both", "Odds", "LogOR", "Concept")
  expect_true(all(res1[,cols] == res2[,cols]))
})