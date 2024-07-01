library(testthat)
library(propr)
library(Rcpp)

test_that("check if odds ratio is computed correctly", {

  # Create two vectors of 1 and 0
  Astar <- c(1, 0, 1, 1, 0, 1, 0, 0, 1, 0)
  Gstar <- c(1, 0, 0, 1, 1, 0, 0, 1, 0, 1)
  
  # Expected values
  a <- 2  # not in A not in G
  b <- 3  # not in A but in G
  c <- 3  # in A but not in G
  d <- 2  # in A in G
  expected_odds_ratio <- (a * d) / (b * c)

  # compute odds ratio table
  or_table <- propr:::getOR(Astar, Gstar)

  # check
  expect_equal(a, or_table[1,1])
  expect_equal(b, or_table[1,2])
  expect_equal(c, or_table[1,3])
  expect_equal(d, or_table[1,4])
  expect_equal(expected_odds_ratio, or_table[1,5])
})

test_that("check if runGraflex produce the expected results", {

  # Create a matrix of 0 and 1
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
                1, 1
              ), nrow = 5, byrow = TRUE)
  colnames(K) <- c("C1", "C2")
  
  # Expected values for C1
  a1 <- 2  # not in A not in G
  b1 <- 2  # not in A but in G
  c1 <- 5  # in A but not in G
  d1 <- 1  # in A in G
  expected_odds_ratio1 <- (a1 * d1) / (b1 * c1)

  # Expected values for C2
  a2 <- 3  # not in A not in G
  b2 <- 1  # not in A but in G
  c2 <- 4  # in A but not in G
  d2 <- 2  # in A in G
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
  
    # Create a matrix of 0 and 1
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
                  1, 1
                ), nrow = 5, byrow = TRUE)
    colnames(K) <- c("C1", "C2")
    
    # compute graflex
    set.seed(0)
    res1 <- propr:::runGraflex(A, K, p=100)
    set.seed(0)
    res2 <- propr:::runGraflex(A, K, p=100)
    set.seed(123)
    res3 <- propr:::runGraflex(A, K, p=100)
  
    # check
    expect_equal(res1, res2)
    expect_false(isTRUE(all.equal(res1, res3)))
})
