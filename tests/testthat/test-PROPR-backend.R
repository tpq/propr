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

test_that("logratio_without_alpha performs log-ratio transformation correctly", {
  # Example data
  example_data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)

  # Apply the logratio_without_alpha function to the example data
  transformed_data <- logratio_without_alpha(example_data, use = 2)

  # Compare the output with the expected output
  expect_identical(transformed_data[1,1], log(1/2))
})
