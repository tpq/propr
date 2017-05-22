if(requireNamespace("SDMTools", quietly = TRUE)){

  x <- abs(rnorm(20))
  w <- abs(rnorm(20))

  test_that("weighted mean matches SDMTools", {

    expect_equal(
      SDMTools::wt.mean(x, w),
      propr:::wtmRcpp(x, w)
    )
  })

  test_that("weighted variance matches SDMTools", {

    expect_equal(
      SDMTools::wt.var(x, w),
      propr:::wtvRcpp(x, w)
    )
  })
}

data(mail)

test_that("lrm without weights matches lrmRcpp", {

  expect_equal(
    propr:::lrmRcpp(mail[]),
    propr:::lrm(mail, mail)
  )
})

test_that("lrv without weights matches vlrRcpp", {

  expect_equal(
    propr:::lltRcpp(propr:::vlrRcpp(mail[])),
    propr:::lrv(mail, mail)
  )
})

if(requireNamespace("limma", quietly = TRUE)){

  data(iris)
  keep <- iris$Species %in% c("setosa", "versicolor")
  counts <- iris[keep, 1:4] * 10
  group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

  test_that("calculateTheta matches calculateThetaW_old", {

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$theta,
      propr::calculateTheta(counts, group, weighted = TRUE)$theta
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm1,
      propr::calculateTheta(counts, group, weighted = TRUE)$lrm1
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm2,
      propr::calculateTheta(counts, group, weighted = TRUE)$lrm2
    )

    expect_equal(
      propr::propd(counts, group, weighted = TRUE)@theta$theta,
      propr::calculateTheta(counts, group, weighted = TRUE)$theta
    )
  })
}
