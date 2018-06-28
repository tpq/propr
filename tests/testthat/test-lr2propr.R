library(propr)

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
      propr:::lr2phi(lr)[1:16],
      phit(ct, symmetrize = FALSE)@matrix[1:16]
    )

    expect_equal(
      propr:::lr2rho(lr)[1:16],
      perb(ct)@matrix[1:16]
    )

    expect_equal(
      propr:::lr2phs(lr)[1:16],
      phis(ct)@matrix[1:16]
    )
  }
})

test_that("lr2 functions handle ivar correctly", {

  expect_equal(
    as.vector(phit(mail, ivar = 2, symmetrize = FALSE)@matrix[, 2]),
    c(Inf, 0, Inf, Inf)
  )

  expect_equal(
    as.vector(perb(mail, ivar = 2)@matrix[, 2]),
    c(0, 1, 0, 0)
  )

  expect_equal(
    as.vector(phis(mail, ivar = 2)@matrix[, 2]),
    c(Inf, 0, Inf, Inf)
  )
})

# Calculate LRV using non-alpha propr
data(mail)
lr <- as.matrix(propr(mail)@logratio)
nonalpha.pr <- propr:::lltRcpp(propr:::lr2vlr(lr))

# Calculate LRV using non-alpha propd
pd <- propd(mail, group = c("A", "A", "A", "B", "B"))
nonalpha.pd <- pd@results$lrv

test_that("lr2vlr matches propd implementation (non-alpha)", {

  expect_equal(
    nonalpha.pr,
    nonalpha.pd
  )
})
