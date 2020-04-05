library(propr)
data(iris)
met.rel <- iris[,1:2]
mic.rel <- iris[,3:4]

# Analyze multi-omics via back-end
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
REL <- cbind(clr(met.rel), clr(mic.rel))
pr.r <- propr:::lr2rho(as.matrix(REL))
colnames(pr.r) <- colnames(REL)
rownames(pr.r) <- colnames(REL)

# Analyze multi-omics with wrapper
ivarNA <- propr(REL, ivar = NA)

test_that("ivar = NA works as expected", {

  expect_equal(
    NA,
    ivarNA@ivar
  )

  expect_equal(
    pr.r,
    ivarNA@matrix
  )

  expect_equal(
    ivarNA@counts,
    ivarNA@logratio
  )
})
