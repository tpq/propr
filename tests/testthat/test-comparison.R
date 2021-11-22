# This is the reference implementation to test against.
ref_impl <- function(values, cutoff, comp_fn) {
  sum(comp_fn(values, cutoff))
}

niters <- 100

test_that("count_less_than works", {
  # This is a ad-hoc property test, without needing to pull in another dependency.
  for (i in 1:niters) {
    vals <- runif(10, -1, 1)
    cutoff <- runif(1, -1, 1)

    expect_equal(count_less_than(vals, cutoff), ref_impl(vals, cutoff, `<`))
    expect_equal(count_greater_than(vals, cutoff), ref_impl(vals, cutoff, `>`))
  }
})
