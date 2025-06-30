library(testthat)
library(propr)

test_that("count_greater_than and count_less_than work properly", {

    # define values
    values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    expect_equal(propr:::count_greater_than(values, 5), 5)
    expect_equal(propr:::count_greater_than(values, 8), 2)
    expect_equal(propr:::count_less_than(values, 5),4)
    expect_equal(propr:::count_less_than(values, 8), 7)

    expect_equal(propr:::count_greater_equal_than(values, 5), 6)
    expect_equal(propr:::count_greater_equal_than(values, 8), 3)
    expect_equal(propr:::count_less_equal_than(values, 5), 5)
    expect_equal(propr:::count_less_equal_than(values, 8), 8)
})

test_that("count_greater_than and count_less_than work properly GPU", {

    # define values
    values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    expect_equal(propr:::count_greater_than(values, 5, use_gpu=TRUE), 5)
    expect_equal(propr:::count_greater_than(values, 8, use_gpu=TRUE), 2)
    expect_equal(propr:::count_less_than(values, 5, use_gpu=TRUE),4)
    expect_equal(propr:::count_less_than(values, 8, use_gpu=TRUE), 7)

    expect_equal(propr:::count_greater_equal_than(values, 5, use_gpu=TRUE), 6)
    expect_equal(propr:::count_greater_equal_than(values, 8, use_gpu=TRUE), 3)
    expect_equal(propr:::count_less_equal_than(values, 5, use_gpu=TRUE), 5)
    expect_equal(propr:::count_less_equal_than(values, 8, use_gpu=TRUE), 8)
})