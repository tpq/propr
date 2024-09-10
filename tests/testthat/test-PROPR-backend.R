library(testthat)
library(propr)

# define data
data <- matrix(c(1:12), ncol=4)
colnames(data) <- c("A", "B", "C", "D")


test_that("index_reference works properly", {

  expect_identical(
    as.integer(index_reference(data, 2)),
    as.integer(2)
  )
  expect_identical(
    as.integer(index_reference(data, c("A", "C"))),
    as.integer(c(1, 3))
  )
  expect_identical(
    as.integer(index_reference(data, "clr")),
    as.integer(c(1, 2, 3, 4))
  )

})

test_that("logratio_without_alpha performs log-ratio transformation correctly", {

  # when ivar is 2
  transformed_data <- logratio_without_alpha(data, use = 2)
  expect_identical(
    as.numeric(round(transformed_data[1,1], 6)), 
    round(log(1/4), 6)
  )
  expect_identical(
    as.numeric(round(transformed_data[2,1], 6)),
    round(log(2/5), 6)
  )

  # when ivar is 1,3
  transformed_data <- logratio_without_alpha(data, use = c(1,3))
  expected <- t(apply(data, 1, function(x) log(x) - mean(log(x[c(1,3)]))))
  expect_true(
    all(round(transformed_data, 6) == round(expected, 6))
  )

  # when ivar is clr
  transformed_data <- logratio_without_alpha(data, use = c(1,2,3,4))
  expected <- t(apply(data, 1, function(x) log(x) - mean(log(x))))
  expect_true(
    all(round(transformed_data, 6) == round(expected, 6))
  )

})
