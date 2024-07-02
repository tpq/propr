library(testthat)
library(propr)

test_that('get propd matrix work',{

  # compute differential proportionality
  x <- iris[,1:4]  # data matrix with 4 variables
  y <- iris[,5]    # group vector
  v <- vegan::rda(log(x[,1]/x[,2]) ~ y)
  pd <- propd(x, as.character(y))
  mat <- getMatrix(pd)

  # check matrix is square and have 4 columns
  expect_equal(ncol(mat), nrow(mat))
  expect_equal(ncol(mat), 4)

  # check diagonal is zero
  expect_equal(mat[1,1], 0)
  expect_equal(mat[2,2], 0)

  # check matrix values
  results <- pd@results
  expect_equal(
    mat[1,2], 
    results[which(results$Partner==2 & results$Pair==1),3]
  )
  expect_equal(
    mat[1,3], 
    results[which(results$Partner==3 & results$Pair==1),3]
  )
})