library(propr)
library(cccrm)
context("linRcpp")

data(mail)
rho <- perb(mail, 3)
swept <- sweep(rho@logratio, 2, apply(rho@logratio, 2, mean), "-")
N <- nrow(swept)

dat <- data.frame("values" = c(swept[, 1], swept[, 2]),
                  "method" = c(rep(1, N), rep(2, N)))
ccc <- cccUst(dat, "values", "method")
names(ccc) <- NULL

lin <- linRcpp(rho@matrix, rho@logratio)

test_that("linRcpp performs Z-transformation correctly (alr)", {

  expect_equal(
    ccc[5],
    lin[2, 1]
  )
})
