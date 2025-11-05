library(testthat)
library(propr)

test_that("3-group RDA agrees with theta", {

  # RDA with 3 groups
  x <- iris[,1:4]
  y <- iris[,5]
  v <- vegan::rda(log(x[,1]/x[,2]) ~ y)
  pd <- propd(x, as.character(y))

  expect_equal(
    1 - sum(v$CCA$eig) / v$tot.chi,
    pd@results[1,"theta"]
  )
})


test_that("2-group RDA agrees with theta", {

  # RDA with 2 groups
  x <- iris[1:100,1:4]
  y <- iris[1:100,5]
  v <- vegan::rda(log(x[,1]/x[,2]) ~ y)
  pd <- propd(x, as.character(y))

  expect_equal(
    1 - sum(v$CCA$eig) / v$tot.chi,
    pd@results[1,"theta"]
  )
})


# data
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

# calculate propd
pd <- propd(counts, group, p = 5)
pd_w <- propd(counts, group, p = 5, weighted = TRUE)

test_that("active theta_e matches calculation using theta_d", {

  n1 <- 50
  n2 <- 50

  expect_equal(
    setActive(pd, what = "theta_e")@results$theta,
    1 - pd@results$theta + pmin((n1-1) * pd@results$lrv1, (n2-1) * pd@results$lrv2) / ((n1+n2-1) * pd@results$lrv)
  )

  # when weighted
  groups <- lapply(unique(group), function(g) g == group)
  ngrp <- length(unique(group))
  # calculate weights, now according to sample reliability weights from limma
  design <-
    stats::model.matrix(~ . + 0, data = as.data.frame(group))

  logX <- log(pd@counts)
  z.geo <- rowMeans(logX)
  z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
  lz.sr <- t(z.lr + mean(z.geo)) #corresponds to log(z.sr) in updateF function

  #use quality weights from limma:
  aw <- limma::arrayWeights(lz.sr, design) 
  W <- t(sweep(matrix(1, nrow(lz.sr), ncol(lz.sr)), 2, aw, `*`)) #get the correct dimensions
  
  ps <- lapply(groups, function(g) propr:::omega(W[g,]))
  names(ps) <- paste0("p", 1:ngrp)
  p <- propr:::omega(W)
  expect_equal(
    setActive(pd_w, what = "theta_e")@results$theta,
    1 - pd_w@results$theta + pmin(ps[[1]] * pd_w@results$lrv1, ps[[2]] * pd_w@results$lrv2) / (p * pd_w@results$lrv)
  )
})

test_that("active theta_f matches calculation using theta_e", {

  expect_equal(
    setActive(pd, what = "theta_f")@results$theta,
    1 - setActive(pd, what = "theta_e")@results$theta
  )

  expect_equal(
    setActive(pd_w, what = "theta_f")@results$theta,
    1 - setActive(pd_w, what = "theta_e")@results$theta
  )
})

test_that("running propd with shrinkage works", {
  pd <- propd(counts, group, p = 5, shrink = TRUE)
  
  # Check that pd is an S4 object and has a "results" slot
  expect_true(isS4(pd))
  expect_true("results" %in% slotNames(pd))
})