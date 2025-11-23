library(testthat)
library(propr)

test_that('get propd matrix work',{

  # compute differential proportionality
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector
  pd <- propd(x, as.character(y))
  pd <- updateF(pd, moderated=T)

  # get results
  results <- pd@results

  # get different theta matrices
  mat_theta <- getMatrix(pd)
  mat_theta_e <- getMatrix(setActive(pd, what='theta_e'))
  mat_theta_mod <- getMatrix(setActive(pd, what='theta_mod'))

  # check matrix is square and have 4 columns
  expect_equal(ncol(mat_theta), nrow(mat_theta))
  expect_equal(ncol(mat_theta), 4)

  # check diagonal is zero
  expect_equal(mat_theta[1,1], 0)
  expect_equal(mat_theta[2,2], 0)

  # check matrix values are correct
  expect_equal(
    mat_theta[1,2], 
    results[which(results$Partner==2 & results$Pair==1),'theta']
  )
  expect_equal(
    mat_theta[2,3], 
    results[which(results$Partner==3 & results$Pair==2),'theta']
  )
  expect_equal(
    mat_theta_e[1,2], 
    results[which(results$Partner==2 & results$Pair==1),'theta_e']
  )
  expect_equal(
    mat_theta_e[2,3], 
    results[which(results$Partner==3 & results$Pair==2),'theta_e']
  )
  expect_equal(
    mat_theta_mod[1,2], 
    results[which(results$Partner==2 & results$Pair==1),'theta_mod']
  )
  expect_equal(
    mat_theta_mod[2,3], 
    results[which(results$Partner==3 & results$Pair==2),'theta_mod']
  )

  # check that theta are different from theta_e and theta_mod
  expect_false(all(mat_theta == mat_theta_e))
  expect_false(all(mat_theta == mat_theta_mod))
})

test_that('get propr matrix work', {

  # create matrix and compute proportionality
  N <- 100
  a <- seq(from = 5, to = 15, length.out = N)
  b <- a * rnorm(N, mean = 1, sd = 0.1)
  c <- rnorm(N, mean = 10)
  d <- rnorm(N, mean = 10)
  e <- rep(10, N)
  X <- data.frame(a, b, c, d, e)
  pr <- propr(X, metric = "rho")

  # check
  expect_equal(
    pr@matrix,
    getMatrix(pr)
  )
})

test_that('results_to_matrix works - when filtered',{

  # compute differential proportionality
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector
  pd <- propd(x, as.character(y))
  pd <- updateF(pd, moderated=T)
  results <- getSignificantResultsFstat(pd, fdr_adjusted=TRUE)
  mat <- results_to_matrix(results)

  # check it is filtered
  expect_true(all(results$FDR <= 0.05))

  # check values are correct
  expect_equal(
    colnames(mat),
    unique(c(results$Pair, results$Partner)),
    info = sprintf(
      "colnames(mat) = c(%s)\nunique(Pair/Partner) = c(%s)",
      paste(colnames(mat), collapse = ", "),
      paste(unique(c(results$Pair, results$Partner)), collapse = ", ")
    )
  )

  theta_msg <- function(i) {
    sprintf(
      "i = %d\nresults$theta[%d] = %s\nmat[%s, %s] = %s",
      i, i,
      results$theta[i],
      results$Pair[i],
      results$Partner[i],
      mat[results$Pair[i], results$Partner[i]]
    )
  }

  expect_equal(
    results$theta[1],
    mat[results$Pair[1], results$Partner[1]],
    info = theta_msg(1)
  )

  expect_equal(
    results$theta[2],
    mat[results$Pair[2], results$Partner[2]],
    info = theta_msg(2)
  )

  expect_equal(
    results$theta[3],
    mat[results$Pair[3], results$Partner[3]],
    info = theta_msg(3)
  )

  # expect_equal(colnames(mat), unique(c(results$Pair, results$Partner)))
  # expect_true(results$theta[1] == mat[results$Pair[1], results$Partner[1]])
  # expect_true(results$theta[2] == mat[results$Pair[2], results$Partner[2]])
  # expect_true(results$theta[3] == mat[results$Pair[3], results$Partner[3]])
})