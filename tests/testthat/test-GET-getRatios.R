library(testthat)
library(propr)

test_that("ratios work when compared with alr", {

  m <- matrix(1:30, 5, 6)
  rat1 <- propr(m, ivar = 1)@logratio[,"6"] # ratio as ALR
  rat2 <- getRatios(propr(m))[,"6/1"] # ratio from getRatios

  expect_equal(
    rat1,
    rat2
  )
})

test_that("ratios function works correctly", {

  # Sample input data
  data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)

  # Test case 1: Without alpha parameter
  expected_output_case1 <- log(data[, 2] / data[, 1])
  result_case1 <- ratios(data)[,1]
  expect_equal(result_case1, expected_output_case1)

  # Test case 2: With alpha parameter
  alpha <- 2
  expected_output_case2 <- (data[, 2]^alpha - data[, 1]^alpha) / alpha
  result_case2 <- ratios(data, alpha = alpha)[,1]
  expect_equal(result_case2, expected_output_case2)
})
