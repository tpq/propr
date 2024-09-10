library(testthat)
library(propr)

test_that("simple_zero_replacement replaces zeros correctly", {

  # Example data with zeros
  data_with_zeros <- matrix(c(1, 2, 3, 0, 5, 6, 0, 8, 9), nrow = 3, byrow = TRUE)

  # Apply the simple_zero_replacement function to the example data
  replaced_data <- simple_zero_replacement(data_with_zeros)

  # Define the expected output after zero replacement
  expected_output <- matrix(c(1, 2, 3, 1, 5, 6, 1, 8, 9), nrow = 3, byrow = TRUE)

  # Compare the output with the expected output
  expect_identical(replaced_data, expected_output)

})