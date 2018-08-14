library(propr)

test_that("ratios function matches ratiosRcpp function", {

  m <- matrix(1:30, 5, 6)
  rat1 <- propr:::ratiosRcpp(m)
  rat2 <- ratios(m)

  expect_equivalent(
    log(rat1),
    rat2
  )

  expect_equivalent(
    ratios(propr(m)@counts)[,"6/2"],
    getRatios(propr(m))[,"6/2"]
  )

  expect_equivalent(
    propr(m, alpha = .5, ivar = 1)@logratio[,"6"],
    getRatios(propr(m, alpha = .5))[,"6/1"]
  )
})
