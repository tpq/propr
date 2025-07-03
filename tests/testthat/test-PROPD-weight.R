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

test_that("test that weights are properly incorporated to lrv", {

    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    v <- limma::voom(t(counts), design = design)
    W <- t(v$weights)

    # calculate lrv using propr
    counts <- as.matrix(counts)
    lrv_w <- propr:::lrv(counts, W, TRUE, NA, counts, W)

    # calculate expected lrv with weights manually
    lrv_e <- c()
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            Wij <- 2 * W[, i] * W[, j] / (W[, i] + W[, j])
            lrv_e <- c(lrv_e, propr:::wtvRcpp(log(counts[, i] / counts[, j]), Wij))
        }
    }

    expect_equal(lrv_w, lrv_e)
    
})

test_that("test that weights are properly incorporated to lrm", {
    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    v <- limma::voom(t(counts), design = design)
    W <- t(v$weights)

    # calculate lrm using propr
    counts <- as.matrix(counts)
    lrm_w <- propr:::lrm(counts, W, TRUE, NA, counts, W)

    # calculate expected lrm with weights manually
    lrm_e <- c()
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            Wij <- 2 * W[, i] * W[, j] / (W[, i] + W[, j])
            lrm_e <- c(lrm_e, propr:::wtmRcpp(log(counts[, i] / counts[, j]), Wij))
        }
    }

    expect_equal(lrm_w, lrm_e)
})
