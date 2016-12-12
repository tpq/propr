library(propr)
context("phis")

data(mail)

test_that("calculating phs from rho matches phis", {

  expect_equal(
    1 - perb(mail)@matrix,
    phis(mail)@matrix
  )
})
