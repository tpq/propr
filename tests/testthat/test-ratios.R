library(propr)

test_that("ratios function matches ratiosRcpp function", {

  m <- matrix(1:30, 5, 6)
  rat1 <- propr:::ratiosRcpp(m)
  rat2 <- ratios(m)

  expect_equivalent(
    rat1,
    rat2
  )
})
