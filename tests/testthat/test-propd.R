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
      propd(counts, group, p = 3)@results$theta
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$theta,
      propd(counts, group, p = 3, weighted = TRUE)@results$theta
    )

    expect_equal(
      propr:::alphaTheta_old(counts, group, alpha = .1)$atheta,
      propd(counts, group, p = 3, alpha = .1)@results$theta
    )

    pdaw <- propd(counts, group, p = 3, weighted = TRUE, alpha = .1)
    expect_equal(
      propr:::alphaThetaW_old(counts, group, alpha = .1, pdaw@weights)$theta,
      pdaw@results$theta
    )
  })

  test_that("propd calculates correct lrm result", {

    expect_equal(
      propr:::calculateTheta_old(counts, group)$lrm1,
      propd(counts, group, p = 3)@results$lrm1
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm1,
      propd(counts, group, p = 3, weighted = TRUE)@results$lrm1
    )

    expect_equal(
      propr:::alphaTheta_old(counts, group, alpha = .1)$alrm1,
      propd(counts, group, p = 3, alpha = .1)@results$lrm1
    )

    pdaw <- propd(counts, group, p = 3, weighted = TRUE, alpha = .1)
    expect_equal(
      propr:::alphaThetaW_old(counts, group, alpha = .1, pdaw@weights)$awlrm1,
      pdaw@results$lrm1
    )
  })
}

test_that("shuffling group labels does not change lrv", {

  expect_equal(
    propr:::calculateTheta(counts[sample(1:100), ], group)$lrv,
    propr:::calculateTheta(counts, group)$lrv
  )

  expect_equal(
    propr:::calculateTheta(counts, group[sample(1:100)])$lrv,
    propr:::calculateTheta(counts, group)$lrv
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
    setActive(pd, what = "theta_e")@results$theta,
    1 - pd@results$theta + pmin((n1-1) * pd@results$lrv1, (n2-1) * pd@results$lrv2) / ((n1+n2-1) * pd@results$lrv)
  )
})

test_that("active theta_f matches calculation using theta_e", {

  expect_equal(
    setActive(pd, what = "theta_f")@results$theta,
    1 - setActive(pd, what = "theta_e")@results$theta
  )
})
