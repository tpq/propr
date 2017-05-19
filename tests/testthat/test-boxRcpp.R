library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

test_that("boxRcpp matches pseudo-lrv from alphaTheta_old", {

  expect_equal(
    propr:::alphaTheta_old(counts, group, .1)$alrv,
    propr:::boxRcpp(as.matrix(counts), .1)
  )
})

test_that("fast alphaTheta matches slow alphaTheta", {

  expect_equal(
    as.vector(propr:::alphaTheta_old(counts, group, .1)$atheta),
    as.vector(propr:::calculateTheta(counts, group, .1)$theta)
  )
})
