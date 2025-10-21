library(testthat)
library(propr)
library(MASS)
library(digest)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

test_that("updateF correctly calculates F-statistics without moderation", {
  
  # Create propd object
  pr <- propd(x, as.character(y), p=10)
  
  # Update F without moderation
  pr_f <- updateF(pr, moderated = FALSE)
  
  # Manual calculation of F-statistic
  N <- length(pr@group)
  K <- length(unique(pr@group))
  expected_fstat <- (N - 2) * (1 - pr@results$theta) / pr@results$theta
  
  # Check F-statistics match
  expect_equal(pr_f@results$Fstat, expected_fstat)
  
  # Check that theta_mod is NA for non-moderated
  expect_true(all(is.na(pr_f@results$theta_mod)))
  
  # Check p-values are calculated
  expect_true("Pval" %in% colnames(pr_f@results))
  expect_true("FDR" %in% colnames(pr_f@results))
  expect_true(all(pr_f@results$Pval >= 0 & pr_f@results$Pval <= 1))
})

test_that("updateF correctly calculates F-statistics with moderation", {
  
  # Create propd object
  pr <- propd(x, as.character(y), p=10)
  
  # Update F with moderation
  pr_f_mod <- updateF(pr, moderated = TRUE)
  
  # Check that F-statistics are calculated
  expect_true("Fstat" %in% colnames(pr_f_mod@results))
  
  # Check that theta_mod is calculated
  expect_true("theta_mod" %in% colnames(pr_f_mod@results))
  expect_true(all(!is.na(pr_f_mod@results$theta_mod)))
  
  # Check p-values are calculated
  expect_true("Pval" %in% colnames(pr_f_mod@results))
  expect_true("FDR" %in% colnames(pr_f_mod@results))
  expect_true(all(pr_f_mod@results$Pval >= 0 & pr_f_mod@results$Pval <= 1))

  # check snapshot values
  expect_equal(
    digest(round(pr_f_mod@results$theta_mod,4)),
    '24f6040126e5fa16fc6c99ffd9aa1959'
  )
})

test_that("updateF correctly calculates F-statistics with moderation and weighted", {
  
  # Create propd object
  pr <- propd(x, as.character(y), p=10, weighted = TRUE)
  
  # Update F with moderation
  pr_f_mod <- updateF(pr, moderated = TRUE)
  
  # Check that F-statistics are calculated
  expect_true("Fstat" %in% colnames(pr_f_mod@results))
  
  # Check that theta_mod is calculated
  expect_true("theta_mod" %in% colnames(pr_f_mod@results))
  expect_true(all(!is.na(pr_f_mod@results$theta_mod)))
  
  # Check p-values are calculated
  expect_true("Pval" %in% colnames(pr_f_mod@results))
  expect_true("FDR" %in% colnames(pr_f_mod@results))
  expect_true(all(pr_f_mod@results$Pval >= 0 & pr_f_mod@results$Pval <= 1))

  # check snapshot values
  expect_equal(
    digest(round(pr_f_mod@results$theta_mod,4)),
    'b89a1343290661cfab2249d506742d1e'
  )
})