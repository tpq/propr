library(propr)
context("lr2")

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
ct <- as.matrix(counts)
lr <- propr:::clrRcpp(ct[])

test_that("lr2 functions match traditional functions", {

  for(i in 1:2){

    expect_equal(
      propr:::lr2vlr(lr),
      propr:::vlrRcpp(ct[])
    )

    expect_equal(
      propr:::lr2phi(lr),
      phit(ct, symmetrize = FALSE)@matrix
    )

    expect_equal(
      propr:::lr2rho(lr),
      perb(ct)@matrix
    )

    expect_equal(
      propr:::lr2phs(lr),
      phis(ct)@matrix
    )
  }
})
