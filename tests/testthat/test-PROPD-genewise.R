library(testthat)
library(propr)

# data
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

# calculate propd and F-statistics (adds FDR column)
pd <- propd(counts, group)
pd <- updateF(pd)

test_that("propdGenewise errors without FDR column", {
  pd_no_fdr <- propd(counts, group)
  expect_error(
    propdGenewise(pd_no_fdr),
    "Please run updateF"
  )
})

test_that("propdGenewise errors with invalid metric", {
  expect_error(
    propdGenewise(pd, metric = "invalid"),
    "'arg' should be one of"
  )
})

test_that("connectivity returns correct structure", {
  res <- propdGenewise(pd, metric = "connectivity")

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("Gene", "connectivity", "FDR"))
  expect_equal(nrow(res), ncol(counts))
  expect_equal(res$Gene, colnames(counts))
})

test_that("wconnectivity returns correct structure", {
  res <- propdGenewise(pd, metric = "wconnectivity")

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("Gene", "wconnectivity", "FDR"))
  expect_equal(nrow(res), ncol(counts))
  expect_equal(res$Gene, colnames(counts))
})

test_that("connectivity counts match manual calculation", {
  fdr_cutoff <- 0.05
  res <- propdGenewise(pd, metric = "connectivity", pairwise_fdr = fdr_cutoff)
  features <- colnames(counts)

  # manually compute connectivity from the results table
  sig <- pd@results[pd@results$FDR > 0 & pd@results$FDR < fdr_cutoff, ]
  manual_conn <- setNames(rep(0L, length(features)), features)
  for (i in seq_len(nrow(sig))) {
    manual_conn[sig$Pair[i]] <- manual_conn[sig$Pair[i]] + 1L
    manual_conn[sig$Partner[i]] <- manual_conn[sig$Partner[i]] + 1L
  }

  expect_equal(res$connectivity, as.integer(manual_conn[features]))
})

test_that("weighted connectivity matches manual calculation", {
  fdr_cutoff <- 0.05
  res <- propdGenewise(pd, metric = "wconnectivity", pairwise_fdr = fdr_cutoff)
  features <- colnames(counts)

  # manually compute weighted connectivity
  sig <- pd@results[pd@results$FDR > 0 & pd@results$FDR < fdr_cutoff, ]
  manual_wconn <- setNames(rep(0, length(features)), features)
  for (i in seq_len(nrow(sig))) {
    w <- 1 / sig$theta[i]
    manual_wconn[sig$Pair[i]] <- manual_wconn[sig$Pair[i]] + w
    manual_wconn[sig$Partner[i]] <- manual_wconn[sig$Partner[i]] + w
  }

  expect_equal(res$wconnectivity, as.numeric(manual_wconn[features]))
})

test_that("strict pairwise_fdr yields zero or fewer connections", {
  res_loose <- propdGenewise(pd, metric = "connectivity", pairwise_fdr = 0.5)
  res_strict <- propdGenewise(pd, metric = "connectivity", pairwise_fdr = 0.001)

  expect_true(all(res_strict$connectivity <= res_loose$connectivity))
})

test_that("pairwise_fdr of 0 gives zero connectivity for all genes", {
  res <- propdGenewise(pd, metric = "connectivity", pairwise_fdr = 0)
  expect_true(all(res$connectivity == 0))

  res_w <- propdGenewise(pd, metric = "wconnectivity", pairwise_fdr = 0)
  expect_true(all(res_w$wconnectivity == 0))
})

test_that("connectivity values are non-negative integers", {
  res <- propdGenewise(pd, metric = "connectivity")
  expect_true(all(res$connectivity >= 0))
  expect_equal(res$connectivity, as.integer(res$connectivity))
})

test_that("wconnectivity values are non-negative", {
  res <- propdGenewise(pd, metric = "wconnectivity")
  expect_true(all(res$wconnectivity >= 0))
})

test_that("default metric is connectivity", {
  res_default <- propdGenewise(pd)
  res_conn <- propdGenewise(pd, metric = "connectivity")
  expect_equal(res_default, res_conn)
})
