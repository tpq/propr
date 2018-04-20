library(propr)

data(mail)

test_that("lrv without weights matches vlrRcpp", {

  expect_equal(
    propr:::lltRcpp(propr:::vlrRcpp(mail[])),
    propr:::lrv(mail, mail)
  )
})

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

if(requireNamespace("limma", quietly = TRUE)){

  test_that("propd returns correct theta result", {

    expect_equal(
      propr:::calculateTheta_old(counts, group)$theta,
      propd(counts, group, p = 3)@theta$theta
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$theta,
      propd(counts, group, p = 3, weighted = TRUE)@theta$theta
    )

    expect_equal(
      propr:::alphaTheta_old(counts, group, alpha = .1)$atheta,
      propd(counts, group, p = 3, alpha = .1)@theta$theta
    )

    pdaw <- propd(counts, group, p = 3, weighted = TRUE, alpha = .1)
    expect_equal(
      propr:::alphaThetaW_old(counts, group, alpha = .1, pdaw@weights)$theta,
      pdaw@theta$theta
    )
  })

  test_that("propd calculates correct lrm result", {

    expect_equal(
      propr:::calculateTheta_old(counts, group)$lrm1,
      propd(counts, group, p = 3)@theta$lrm1
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm1,
      propd(counts, group, p = 3, weighted = TRUE)@theta$lrm1
    )

    expect_equal(
      propr:::alphaTheta_old(counts, group, alpha = .1)$alrm1,
      propd(counts, group, p = 3, alpha = .1)@theta$lrm1
    )

    pdaw <- propd(counts, group, p = 3, weighted = TRUE, alpha = .1)
    expect_equal(
      propr:::alphaThetaW_old(counts, group, alpha = .1, pdaw@weights)$awlrm1,
      pdaw@theta$lrm1
    )
  })
}

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
theta <- propr:::calculateTheta_old(counts, group)
ptheta <- propr:::permuteTheta_prime(counts, group, p = 5)
pt <- propr:::calculateFDR(theta, ptheta, cutoff = seq(.95, 1, .01))

set.seed(1)
pd <- propd(counts, group, p = 5)
pd <- updateCutoffs(pd, cutoff = seq(.95, 1, .01))

test_that("propd FDR mirrors permuteTheta_prime", {

  expect_equal(
    pt$FDR,
    pd@fdr$FDR
  )
})

n1 <- 50
n2 <- 50

test_that("active theta_e matches calculation using theta_d", {

  expect_equal(
    setActive(pd, what = "theta_e")@theta$theta,
    1 - pd@theta$theta + pmin((n1-1) * pd@theta$lrv1, (n2-1) * pd@theta$lrv2) / ((n1+n2-1) * pd@theta$lrv)
  )
})

test_that("active theta_f matches calculation using theta_e", {

  expect_equal(
    setActive(pd, what = "theta_f")@theta$theta,
    1 - setActive(pd, what = "theta_e")@theta$theta
  )
})
