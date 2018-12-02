library(propr)
data(caneToad.counts)
data(caneToad.groups)

df <- caneToad.counts[,1:500]
grp <- caneToad.groups
A <- propr(df, "phs")

test_that("getMatrix agrees with original matrix", {

  expect_equal(
    A@matrix,
    getMatrix(A)
  )
})
