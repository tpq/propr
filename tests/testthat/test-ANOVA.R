library(propr)
data(iris)

x <- iris[,1:4]
y <- iris[,5]
v <- vegan::rda(log(x[,1]/x[,2]) ~ y)
pd <- propd(x, as.character(y))

test_that("3-group RDA agrees with theta", {

  expect_equal(
    1 - sum(v$CCA$eig) / v$tot.chi,
    pd@results[1,"theta"]
  )
})

x <- iris[1:100,1:4]
y <- iris[1:100,5]
v <- vegan::rda(log(x[,1]/x[,2]) ~ y)
pd <- propd(x, as.character(y))

test_that("2-group RDA agrees with theta", {

  expect_equal(
    1 - sum(v$CCA$eig) / v$tot.chi,
    pd@results[1,"theta"]
  )
})
