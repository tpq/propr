library(propr)
context("lrmodel")

library(propr)
data(mail)
rho <- perb(mail)

test_that("lrmodel correctly builds and deploys rule", {

  model <- modelCLR(mail)
  lr <- predict(model, mail)

  expect_equivalent(
    model@data,
    lr
  )

  expect_equivalent(
    model@data,
    rho@logratio
  )

  model2 <- modelCLR(t(mail), MARGIN = 2)
  lr2 <- predict(model2, t(mail))

  expect_equivalent(
    model2@data,
    lr2
  )

  expect_equivalent(
    model2@data,
    t(rho@logratio)
  )
})
