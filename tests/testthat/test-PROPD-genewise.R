library(testthat)
library(propr)

# ---- shared test data (simulated, for genewise tests) ----
set.seed(42)

# simulate 100 genes, 40 samples, 2 groups
n_samples <- 40
n_genes   <- 100
group     <- rep(c("A", "B"), each = n_samples / 2)

# base counts with some differential proportionality signal
counts_base <- matrix(rnbinom(n_samples * n_genes, mu = 100, size = 10),
                      nrow = n_samples, ncol = n_genes)
colnames(counts_base) <- paste0("gene", seq_len(n_genes))
rownames(counts_base) <- paste0("sample", seq_len(n_samples))

# add signal: make a few genes differentially proportional in group B
signal_genes <- sample(seq(100))[1:10]
counts_base[group == "B", signal_genes] <- matrix(
  rnbinom(sum(group == "B") * length(signal_genes), mu = 300, size = 5),
  nrow = sum(group == "B"))

propd_testdata <- list(counts = counts_base, group = group)

# ---- Generate propd object ------------

counts <- propd_testdata$counts
group  <- propd_testdata$group

pd        <- propd(counts, group)
pd_with_f <- updateF(pd)

# Run genewise once; reuse across tests to keep the suite fast.
res <- propdGenewise(pd_with_f)

# ---- error conditions ----

test_that("propdGenewise errors without FDR column", {
  expect_error(propdGenewise(pd), "Please run updateF")
})

test_that("propdGenewise errors for alpha != NA", {
  pd_alpha <- propd(counts, group, alpha = 0.5)
  pd_alpha <- updateF(pd_alpha)
  expect_error(propdGenewise(pd_alpha), "only works for alpha = NA")
})

test_that("propdGenewise errors for more than 2 groups", {
  pd3 <- propd(iris[, 1:4], as.character(iris[, 5]))
  pd3 <- updateF(pd3)
  expect_error(propdGenewise(pd3), "only works for 2 groups")
})

# ---- output structure ----

test_that("propdGenewise returns a data frame with correct columns", {
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res),
               c("id", "lfc", "lrmD", "connectivity", "ES", "padj", "ES_batch"))
})

test_that("propdGenewise has one row per gene", {
  expect_equal(nrow(res), ncol(counts))
  expect_equal(res$id, colnames(counts))
})

# ---- connectivity ----

test_that("connectivity matches manual calculation", {
  fdr_cutoff <- 0.05
  res_c <- propdGenewise(pd_with_f, pairwise_fdr = fdr_cutoff)
  features <- colnames(counts)

  sig <- pd_with_f@results[
    pd_with_f@results$FDR > 0 & pd_with_f@results$FDR < fdr_cutoff, ]
  manual_conn <- setNames(rep(0L, length(features)), features)
  for (i in seq_len(nrow(sig))) {
    manual_conn[sig$Pair[i]]    <- manual_conn[sig$Pair[i]]    + 1L
    manual_conn[sig$Partner[i]] <- manual_conn[sig$Partner[i]] + 1L
  }

  expect_equal(res_c$connectivity, as.integer(manual_conn[features]))
})

test_that("connectivity values are non-negative integers", {
  expect_true(all(res$connectivity >= 0))
  expect_equal(res$connectivity, as.integer(res$connectivity))
})

test_that("stricter pairwise_fdr yields fewer or equal connections", {
  res_loose  <- propdGenewise(pd_with_f, pairwise_fdr = 0.5)
  res_strict <- propdGenewise(pd_with_f, pairwise_fdr = 0.001)
  expect_true(all(res_strict$connectivity <= res_loose$connectivity))
})

test_that("pairwise_fdr = 0 gives zero connectivity for all genes", {
  res_zero <- propdGenewise(pd_with_f, pairwise_fdr = 0)
  expect_true(all(res_zero$connectivity == 0))
})

# ---- lfc ----

test_that("lfc is computed correctly via CLR", {
  res <- propdGenewise(pd_with_f)

  # manually compute CLR-based LFC
  ct <- as.matrix(counts)
  if (any(ct == 0)) ct <- ct + 1

  g1 <- group == "A"
  g2 <- group == "B"

  clr1 <- propr:::logratio(ct[g1,], 'clr', NA)
  clr2 <- propr:::logratio(ct[g2,], 'clr', NA)
  manual_lfc <- (colMeans(clr1 - clr2)) / log(2)

  expect_equal(res$lfc, as.numeric(manual_lfc))
})

# ---- lrmD ----

test_that("lrmD is NA when gene has no significant connections", {
  # use very strict cutoff so no pairs are significant
  res <- propdGenewise(pd_with_f, pairwise_fdr = 0)
  expect_true(all(is.na(res$lrmD)))
})

# ---- ES / padj ----

test_that("ES values are finite", {
  expect_true(all(is.finite(res$ES)))
})

test_that("padj values are between 0 and 1", {
  expect_true(all(res$padj >= 0 & res$padj <= 1, na.rm = TRUE))
})

test_that("ES_batch values are finite", {
  expect_true(all(is.finite(res$ES_batch)))
})

test_that("propdGenewise is reproducible with same inputs", {
  res1 <- propdGenewise(pd_with_f)
  res2 <- propdGenewise(pd_with_f)
  expect_equal(res1$padj, res2$padj)
  expect_equal(res1$ES_batch, res2$ES_batch)
})

test_that("ES values are non-negative when scoreType is pos", {
  expect_true(all(res$ES >= 0))
})

# ---- .compute_es internal function ----

test_that(".compute_es returns expected structure", {
  es <- propr:::.compute_es(pd_with_f@results, colnames(counts))

  expect_s3_class(es, "data.frame")
  expect_equal(colnames(es), c("es", "es_pos", "es_neg"))
  expect_equal(nrow(es), ncol(counts))
  expect_true(all(is.finite(es$es)))
  expect_true(all(is.finite(es$es_pos)))
  expect_true(all(is.finite(es$es_neg)))
})

test_that(".compute_es: es_pos is non-negative and es_neg is non-positive", {
  es <- propr:::.compute_es(pd_with_f@results, colnames(counts))
  expect_true(all(es$es_pos >= 0))
  expect_true(all(es$es_neg <= 0))
})

test_that(".compute_es: es is the larger-magnitude of es_pos and es_neg", {
  es <- propr:::.compute_es(pd_with_f@results, colnames(counts))
  expected_es <- ifelse(abs(es$es_pos) >= abs(es$es_neg), es$es_pos, es$es_neg)
  expect_equal(es$es, expected_es)
})
