library(testthat)
library(propr)

# define data matrix
set.seed(123)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

# old code
updateCutoffs_old <- function(object, cutoff, ncores=1){

  getFdrRandcounts <- function(ct.k){

    pr.k <- suppressMessages(
      propr::propr(ct.k, object@metric, ivar = object@ivar, alpha = object@alpha, p = 0))

    # Vector of propr scores for each pair of taxa.
    pkt <- pr.k@results$propr

    # Find number of permuted theta less than cutoff
    sapply(FDR$cutoff, function(cut){
      if(object@metric == "rho" | object@metric == "cor"){
        propr:::count_greater_than(pkt, cut)
      }else{ # phi & phs
        propr:::count_less_than(pkt, cut)
      }
    })
  }

  if(object@metric == "rho"){

    message("Alert: Estimating FDR for largely positive proportional pairs only.")
  }

  if(object@metric == "phi"){

    warning("We recommend using the symmetric phi 'phs' for FDR permutation.")
  }

  if(identical(object@permutes, list(NULL))) stop("Permutation testing is disabled.")

  # Let NA cutoff skip function
  if(identical(cutoff, NA)) return(object)

  # Set up FDR cutoff table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
  FDR$cutoff <- cutoff
  p <- length(object@permutes)

  if(ncores > 1){
    skip(message="WARN:Skipping cannot run on local CI")
    packageCheck("parallel")

    # Set up the cluster and require propr
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterEvalQ(cl, requireNamespace(propr, quietly = TRUE))

    # Each element of this list will be a vector whose elements
    # are the count of theta values less than the cutoff.
    randcounts <- parallel::parLapply(
      cl = cl,
      X = object@permutes,
      fun = getFdrRandcounts
    )

    # Sum across cutoff values
    FDR$randcounts <- apply(as.data.frame(randcounts), 1, sum)

    # Explicitly stop the cluster.
    parallel::stopCluster(cl)

  }else{

    # Calculate propr for each permutation -- NOTE: `select` and `subset` disable permutation testing
    for(k in 1:p){

      # Calculate propr exactly based on @metric, @ivar, and @alpha
      ct.k <- object@permutes[[k]]
      pr.k <- suppressMessages(
        propr(ct.k, object@metric, ivar = object@ivar, alpha = object@alpha, p = 0))
      pkt <- pr.k@results$propr

      # Find number of permuted theta less than cutoff
      for(cut in 1:nrow(FDR)){ # randcounts as cumsum

        # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
        if(object@metric == "rho" | object@metric == "cor"){
          FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + propr:::count_greater_than(pkt, FDR[cut, "cutoff"])
        }else{ # phi & phs
          FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + propr:::count_less_than(pkt, FDR[cut, "cutoff"])
        }
      }
    }
  }

  # Calculate FDR based on real and permuted tallys
  FDR$randcounts <- FDR$randcounts / p # randcounts as mean

  for(cut in 1:nrow(FDR)){

    # Count positives as rho > cutoff, cor > cutoff, phi < cutoff, phs < cutoff
    if(object@metric == "rho" | object@metric == "cor"){
      FDR[cut, "truecounts"] <- sum(object@results$propr > FDR[cut, "cutoff"])
    }else{ # phi & phs
      FDR[cut, "truecounts"] <- sum(object@results$propr < FDR[cut, "cutoff"])
    }

    FDR[cut, "FDR"] <- FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
  }

  # Initialize @fdr
  object@fdr <- FDR

  return(object)
}

test_that("updateCutoffs_old and updateCutoffs.propr agree - rho", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "rho", p=10, permutation_option = "sample-wise")
  set.seed(0)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)
  set.seed(0)
  pr_old <- updateCutoffs_old(pr, cutoff = pr@fdr$cutoff)

  # check that the two methods agree
  expect_equal(pr@fdr, pr_old@fdr)
})

test_that("updateCutoffs_old and updateCutoffs.propr agree - phi", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "phi", p=10, permutation_option = "sample-wise")
  set.seed(0)
  pr <- updateCutoffs(pr, number_of_cutoffs=10)
  set.seed(0)
  pr_old <- updateCutoffs_old(pr, cutoff = pr@fdr$cutoff)

  # check that the two methods agree
  expect_equal(pr@fdr, pr_old@fdr)
})

test_that("updateCutoffs_old and updateCutoffs.propr agree - cor", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "cor", p=10, permutation_option = "sample-wise")
  set.seed(0)
  pr <- updateCutoffs(pr, number_of_cutoffs=10, tails="right")
  set.seed(0)
  pr_old <- updateCutoffs_old(pr, cutoff = pr@fdr$cutoff)

  # check that the two methods agree
  expect_equal(pr@fdr, pr_old@fdr)
})

test_that("updateCutoffs.propr properly set up cutoffs", {

  # get propr object
  pr <- propr(X, metric = "pcor.bshrink", p=10)

  # get cutoffs
  values <- pr@results$propr
  cutoffs_right <- as.numeric( quantile(values[values >= 0], probs = seq(0, 1, length.out = 10)) )
  cutoffs_both <- as.numeric( quantile(abs(values), probs = seq(0, 1, length.out = 10)) )

  # check that cutoffs are properly defined
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10, tails="right")@fdr$cutoff,
    cutoffs_right
  )
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10, tails="both")@fdr$cutoff,
    cutoffs_both
  )
  expect_equal(
    updateCutoffs(pr, number_of_cutoffs=10)@fdr$cutoff,
    cutoffs_right
  )
  expect_equal(
    updateCutoffs(pr, custom_cutoffs=cutoffs_right)@fdr,
    updateCutoffs(pr, number_of_cutoffs=10)@fdr
  )
  expect_equal(
    updateCutoffs(pr, custom_cutoffs=cutoffs_both, tails='both')@fdr,
    updateCutoffs(pr, number_of_cutoffs=10, tails='both')@fdr
  )
})

test_that("updateCutoffs.propr properly calculates truecounts", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr_right <- updateCutoffs(pr, number_of_cutoffs=10, tails='right')
  pr_both <- updateCutoffs(pr, number_of_cutoffs=10, tails='both')
  pr_default <- updateCutoffs(pr, number_of_cutoffs=10)

  # get truecounts
  truecounts_right <- sapply(
    pr@fdr$cutoff, 
    function(cut) sum(pr@results$propr > cut)
  )
  truecounts_both <- sapply(
    pr@fdr$cutoff, 
    function(cut) sum(abs(pr@results$propr) > cut)
  )

  # check that truecounts are properly defined
  expect_true(all(pr_right@fdr$truecounts == truecounts_right))
  expect_true(all(pr_both@fdr$truecounts == truecounts_both))
  expect_true(all(pr_default@fdr$truecounts == truecounts_both))
})

test_that("updateCutoffs.propr properly calculates randcounts", {

  # get propr object and update cutoffs
  pr <- propr(X, metric = "pcor.bshrink", p=10)
  pr <- updateCutoffs(pr, number_of_cutoffs=10, tails='both')

  # get permuted values
  randcounts <- rep(0, 10)
  for (k in 1:10){
    ct.k <- pr@permutes[[k]]
    pr.k <- suppressMessages(propr(
      ct.k,
      pr@metric,
      ivar = pr@ivar,
      alpha = pr@alpha,
      p = 0
    ))
    pkt <- pr.k@results$propr
    pkt <- abs(pkt)
    for (cut in 1:length(pr@fdr$cutoff)){
      cutoff <- pr@fdr$cutoff[cut]
      randcounts[cut] <- randcounts[cut] + sum(pkt > cutoff)
    }
  }
  randcounts <- randcounts / 10

  # check that the permutation tests work properly
  expect_equal(pr@fdr$randcounts, randcounts)
})

test_that("updateCutoffs.propr is reproducible when seed is set", {
  
    # get propr object and update cutoffs
    set.seed(0)
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    set.seed(0)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
    pr3 <- propr(X, metric = "pcor.bshrink", p=10)
    pr3 <- updateCutoffs(pr3, number_of_cutoffs=10)
  
    # check that fdr are the same only when seed is the same
    expect_equal(pr1@fdr, pr2@fdr)
    expect_false(isTRUE(all.equal(pr1@fdr, pr3@fdr)))
})

test_that("updateCutoffs.propr works when ncores > 1", {
    skip(message="WARN:Skipping cannot run on local CI")
    # get propr object and update cutoffs
    set.seed(0)
    pr1 <- propr(X, metric = "pcor.bshrink", p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10, ncores=1)
    set.seed(0)
    pr2 <- propr(X, metric = "pcor.bshrink", p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10, ncores=2)

    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})