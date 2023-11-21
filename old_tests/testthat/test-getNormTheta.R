library(propr)

test_that("getNormTheta works as expected", {

  x <- matrix(1:30, 5, 6)
  y <- c("A", "A", "A", "B", "B")
  pd <- propd(x, y)

  # Get per-feature theta for first feature manually
  thetaFor1 <- pd@results[pd@results$Pair == 1, ]
  theta.1 <- thetaFor1$theta
  names(theta.1) <- colnames(pd@counts)[thetaFor1$Partner]

  # Get per-feature theta for first feature
  norm.1 <- getNormTheta(pd, x[,1])

  expect_equal(
    theta.1,
    norm.1[-1] # first entry is DP with self = 1
  )
})
