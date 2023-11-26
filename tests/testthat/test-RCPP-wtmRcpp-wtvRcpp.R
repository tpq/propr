# library(testthat)
# library(propr)
#
# if(requireNamespace("SDMTools", quietly = TRUE)){
#
#   x <- abs(rnorm(20))
#   w <- abs(rnorm(20))
#
#   test_that("weighted mean matches SDMTools", {
#
#     expect_equal(
#       SDMTools::wt.mean(x, w),
#       propr:::wtmRcpp(x, w)
#     )
#   })
#
#   test_that("weighted variance matches SDMTools", {
#
#     expect_equal(
#       SDMTools::wt.var(x, w),
#       propr:::wtvRcpp(x, w)
#     )
#   })
# }
