library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

test_that("propd returns correct thetea result", {

  expect_equal(
    propriety:::calculateTheta(counts, group)$theta,
    propd(counts, group, p = 3)@theta$theta
  )

  expect_equal(
    propriety:::alphaTheta(counts, group, alpha = .1)$theta,
    propd(counts, group, p = 3, alpha = .1)@theta$theta
  )
})

test_that("shuffling group labels does not change lrv", {

  expect_equal(
    calculateTheta(counts[sample(1:100), ], group)$lrv,
    calculateTheta(counts, group)$lrv
  )

  expect_equal(
    calculateTheta(counts, group[sample(1:100)])$lrv,
    calculateTheta(counts, group)$lrv
  )
})

set.seed(1)
theta <- propriety:::calculateTheta(counts, group)
ptheta <- propriety:::permuteTheta_prime(counts, group, p = 5)
pt <- propriety:::calculateFDR(theta, ptheta, cutoff = seq(.95, 1, .01))

set.seed(1)
pd <- propd(counts, group, p = 5, cutoff = seq(.95, 1, .01))

test_that("propd FDR mirrors permuteTheta_prime", {

  expect_equal(
    pt$FDR,
    pd@fdr$FDR
  )
})
