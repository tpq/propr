library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

rho <- perb(counts)@matrix
diag(rho) <- 0

test_that("half-matrix correctly turned into matrix", {

  expect_equal(
    rho,
    propr:::half2mat(propr:::lltRcpp(rho))
  )
})
