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

# Test updateCutoffs for the ivar = NA wrapper
set.seed(1)
pr_auto <- propr(iris[,1:4])
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
myCLR <- clr(iris[,1:4])
set.seed(1)
pr_manual <- propr(myCLR, ivar = NA)
pr_auto <- updateCutoffs(pr_auto, cutoff = seq(0, 1, .05))
pr_manual <- updateCutoffs(pr_manual, cutoff = seq(0, 1, .05))

test_that("ivar = NA will work with FDR", {

  expect_equal(
    pr_auto@fdr,
    pr_manual@fdr
  )
})
