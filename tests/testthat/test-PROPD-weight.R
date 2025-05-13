library(testthat)
library(propr)

# data
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

test_that("test that providing weights to propd works", {

    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    v <- limma::voom(t(counts), design = design)
    W <- t(v$weights)

    # calculate propr 
    pd_w <- propd(counts, group, weighted = TRUE)
    pd_w2 <- propd(counts, group, weighted = TRUE, weights = W)

    expect_equal(pd_w@weights, pd_w2@weights)
    expect_equal(pd_w@results, pd_w2@results)
})