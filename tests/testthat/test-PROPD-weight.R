library(testthat)
library(propr)

fetch_weights <- function(counts, design){
    #this function calculates the new default weights, i.e., sample reliability weights using limma
    logX <- log(counts)
    z.geo <- rowMeans(logX)
    z.lr <- as.matrix(sweep(logX, 1, z.geo, "-"))
    lz.sr <- t(z.lr + mean(z.geo)) #corresponds to log(z.sr) in updateF function

    #use quality weights from limma:
    aw <- limma::arrayWeights(lz.sr, design) 
    W <- t(sweep(matrix(1, nrow(lz.sr), ncol(lz.sr)), 2, aw, `*`)) #get the correct dimensions
    return(W)
}  

# data
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")


test_that("test that weights are properly incorporated to lrv", {

    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    W <- fetch_weights(counts,design)

    # calculate lrv using propr
    counts <- as.matrix(counts)
    lrv_w <- propr:::lrv(counts, W, TRUE, NA, counts, W)

    # calculate expected lrv with weights manually
    lrv_e <- c()
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            Wij <- 2 * W[, i] * W[, j] / (W[, i] + W[, j])
            x <- log(counts[, i] / counts[, j])
            xbar <- sum(x * Wij) / sum(Wij)
            lrv <- sum(Wij * (x - xbar)^2) / (sum(Wij) - sum(Wij^2) / sum(Wij))
            lrv_e <- c(lrv_e, lrv)
            #lrv_e <- c(lrv_e, propr:::wtvRcpp(log(counts[, i] / counts[, j]), Wij))  # this is the same as above
        }
    }

    expect_equal(lrv_w, lrv_e)
    
})

test_that("test that weights are properly incorporated to lrm", {
    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    W <- fetch_weights(counts,design)

    # calculate lrm using propr
    counts <- as.matrix(counts)
    lrm_w <- propr:::lrm(counts, W, TRUE, NA, counts, W)

    # calculate expected lrm with weights manually
    lrm_e <- c()
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            Wij <- 2 * W[, i] * W[, j] / (W[, i] + W[, j])
            x <- log(counts[, i] / counts[, j])
            xbar <- sum(x * Wij) / sum(Wij)
            lrm_e <- c(lrm_e, xbar)
        }
    }

    expect_equal(lrm_w, lrm_e)
})

test_that("test that weights are properly incorporated to omega" ,{
    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    W <- fetch_weights(counts,design)

    # calculate omega using propr
    counts <- as.matrix(counts)
    omega_w <- propr:::omega(W)

    # calculate expected omega with weights manually
    omega_e <- c()
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            Wij <- 2 * W[, i] * W[, j] / (W[, i] + W[, j])
            omega_e <- c(omega_e, sum(Wij) - sum(Wij^2) / sum(Wij))
        }
    }

    expect_equal(omega_w, omega_e)
})

test_that("test that weights are properly incorporated to theta", {
    # get weights
    design <- stats::model.matrix(~ . + 0, data = as.data.frame(group))
    W <- fetch_weights(counts,design)

    # calculate theta using propr
    counts <- as.matrix(counts)
    theta_w <- propr:::calculate_theta(counts, group, weighted=TRUE)$theta

    # calculate expected theta with weights manually
    theta_e <- c()
    groups <- lapply(unique(group), function(g) g == group)
    Wfull <- W
    W1 <- W[groups[[1]],]
    W2 <- W[groups[[2]],]
    counts1 <- counts[groups[[1]],]
    counts2 <- counts[groups[[2]],]
    for (i in 2:ncol(counts)) {
        for (j in 1:(i-1)) {
            # calculate lrv and omega for group 1
            Wij1 <- 2 * W1[, i] * W1[, j] / (W1[, i] + W1[, j])
            lrv1 <- propr:::wtvRcpp(log(counts1[, i] / counts1[, j]), Wij1)
            omega1 <- sum(Wij1) - sum(Wij1^2) / sum(Wij1)

            # calculate lrv and omega for group 2
            Wij2 <- 2 * W2[, i] * W2[, j] / (W2[, i] + W2[, j])
            lrv2 <- propr:::wtvRcpp(log(counts2[, i] / counts2[, j]), Wij2)
            omega2 <- sum(Wij2) - sum(Wij2^2) / sum(Wij2)

            # calculate lrv and omega between groups
            Wij <- 2 * Wfull[, i] * Wfull[, j] / (Wfull[, i] + Wfull[, j])
            lrv <- propr:::wtvRcpp(log(counts[, i] / counts[, j]), Wij)
            omega <- sum(Wij) - sum(Wij^2) / sum(Wij)

            # calculate theta
            theta <- (omega1 * lrv1 + omega2 * lrv2) / (omega * lrv)
            bad <- is.na(lrv1) || is.na(lrv2) || is.na(lrv) || (lrv == 0)
            if (bad || !is.finite(theta)) {
                theta <- 1
            }
            theta_e <- c(theta_e, theta)
        }
    }

    expect_equal(theta_w, theta_e)
})
