library(propr)
context("phis")

data(mail)
lr <- propr:::clrRcpp(mail[])

test_that("calculating phs from rho matches phis", {

  expect_equal(
    var(lr[, 1] - lr[, 2]) / var(lr[, 1] + lr[, 2]),
    phis(mail)@matrix[1, 2]
  )

  expect_equal(
    propr:::lr2phs(lr),
    (1 - propr:::lr2rho(lr)) / (1 + propr:::lr2rho(lr))
  )
})
