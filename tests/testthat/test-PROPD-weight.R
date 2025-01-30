library(testthat)
library(propr)

# data
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

test_that("test that providing weights to propd works", {
    pd_w <- propd(counts, group, p = 5, weighted = TRUE)
    w <- pd_w@weights
    pd_w2 <- propd(counts, group, p = 5, weighted = TRUE, weights = w)

    expect_equal(pd_w@weights, pd_w2@weights)
    expect_equal(pd_w@results, pd_w2@results)
})