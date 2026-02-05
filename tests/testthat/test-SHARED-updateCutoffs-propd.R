library(testthat)
library(propr)
library(MASS)

# define data
data(crabs)
x <- crabs[,4:8]  # data matrix with 5 variables
y <- crabs[,1]    # group vector

# old code
updateCutoffs_old <- function(object, cutoff = seq(.05, .95, .3)){

  if(identical(object@permutes, data.frame())) stop("Permutation testing is disabled.")

  # Let NA cutoff skip function
  if(identical(cutoff, NA)) return(object)

  # Set up FDR cutoff table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", "FDR")
  FDR$cutoff <- cutoff
  p <- ncol(object@permutes)
  lrv <- object@results$lrv

  # Use calculateTheta to permute active theta
  for(k in 1:p){

    # Tally k-th thetas that fall below each cutoff
    shuffle <- object@permutes[, k]

    if(object@active == "theta_mod"){

      # Calculate theta_mod with updateF (using i-th permuted object)
      if(is.na(object@Fivar)) stop("Please re-run 'updateF' with 'moderation = TRUE'.")
      propdi <- suppressMessages(
        propd(object@counts[shuffle, ], group = object@group, alpha = object@alpha, p = 0,
              weighted = object@weighted))
      propdi <- suppressMessages(
        updateF(propdi, moderated = TRUE, ivar = object@Fivar))
      pkt <- propdi@results$theta_mod

    }else{

      # Calculate all other thetas directly (using calculateTheta)
      pkt <- suppressMessages(
        calculate_theta(object@counts[shuffle, ], object@group, object@alpha, lrv,
                       only = object@active, weighted = object@weighted))
    }

    # Find number of permuted theta less than cutoff
    for(cut in 1:nrow(FDR)){ # randcounts as cumsum
      FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + sum(pkt < FDR[cut, "cutoff"])
    }
  }

  # Calculate FDR based on real and permuted tallys
  FDR$randcounts <- FDR$randcounts / p # randcounts as mean
  for(cut in 1:nrow(FDR)){
    FDR[cut, "truecounts"] <- sum(object@results$theta < FDR[cut, "cutoff"])
    FDR[cut, "FDR"] <- FDR[cut, "randcounts"] / FDR[cut, "truecounts"]
  }

  # Initialize @fdr
  object@fdr <- FDR

  return(object)
}

test_that("updateCutoffs_old and updateCutoffs.propd agree", {

  # get propd object and update cutoffs
  pd <- propd(x, as.character(y), p=10)
  set.seed(0)
  pd <- updateCutoffs(pd, number_of_cutoffs=10)
  set.seed(0)
  pd_old <- updateCutoffs_old(pd, cutoff = pd@fdr$cutoff)

  # check that the two methods agree
  expect_equal(pd@fdr, pd_old@fdr)
})

test_that("updateCutoffs_old and updateCutoffs.propd agree - when using theta_mod", {

  # use this data instead, otherwise the limma voom will return Inf prior dfs, probably because crab data is not count data
  x2 <- iris[,1:4]
  y2 <- iris[,5]

  # get propd object and update cutoffs
  pd <- propd(x2, as.character(y2), p=10)
  pd <- updateF(pd, moderated = TRUE)
  pd <- setActive(pd, "theta_mod")
  set.seed(0)
  pd <- updateCutoffs(pd, number_of_cutoffs=10)
  set.seed(0)
  pd_old <- updateCutoffs_old(pd, cutoff = pd@fdr$cutoff)

  # check that the two methods agree
  expect_equal(pd@fdr, pd_old@fdr)
})

test_that("updateCutoffs.propd properly set up cutoffs", {

  # get propd object and update cutoffs
  pd <- propd(x, as.character(y), p=10)
  pd1 <- updateCutoffs(pd, number_of_cutoffs=10)

  # get cutoffs
  cutoffs <- as.numeric( quantile(pd@results$theta, probs = seq(0, 1, length.out = 10)) )

  # run with cutoffs
  pd2 <- updateCutoffs(pd, custom_cutoffs=cutoffs)

  # check that cutoffs are properly defined
  expect_equal(pd1@fdr$cutoff, cutoffs)

  # check that the two calls agree
  expect_equal(pd1@fdr, pd2@fdr)
})

test_that("updateCutoffs.propd properly calculates truecounts", {

  # get propd object and update cutoffs
  pd <- propd(x, as.character(y), p=10)
  pd <- updateCutoffs(pd, number_of_cutoffs=10)

  # get truecounts
  truecounts <- sapply(
    pd@fdr$cutoff, 
    function(cut) sum(pd@results$theta < cut)
  )
  truecounts_manual <- c(0:9)

  # check that truecounts are properly defined
  expect_equal(pd@fdr$truecounts, truecounts)
  expect_equal(pd@fdr$truecounts, truecounts_manual)
})

test_that("updateCutoffs.propd properly calculates randcounts", {

  # get propd object and update cutoffs
  set.seed(0)
  pd <- propd(x, as.character(y), p=10)
  pd <- updateCutoffs(pd, number_of_cutoffs=10)

  # get permuted values
  randcounts <- rep(0, 10)
  for (k in 1:10){
    shuffle <- pd@permutes[,k]
    ct.k <- pd@counts[shuffle,]
    pd.k <- suppressMessages(propd(
      ct.k,
      group = pd@group,
      alpha = pd@alpha,
      p = 0,
      weighted = pd@weighted
    ))
    pkt <- pd.k@results$theta
    for (cut in 1:length(pd@fdr$cutoff)){
        randcounts[cut] <- randcounts[cut] + sum(pkt < pd@fdr$cutoff[cut])
    }
  }
  randcounts <- randcounts / 10
  randcounts_manual <- c(0, 0, 0, 0, 0, 0, 0, 0.1, 0.1, 5.4)

  # check that the permutation tests work properly
  expect_equal(pd@fdr$randcounts, randcounts)
  expect_equal(round(pd@fdr$randcounts, 1), randcounts_manual)
})

test_that("updateCutoffs.propd is reproducible when seed is set", {
  
    # get propd object and update cutoffs
    set.seed(0)
    pr1 <- propd(x, as.character(y), p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10)
    set.seed(0)
    pr2 <- propd(x, as.character(y), p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10)
    pr3 <- propd(x, as.character(y), p=10)
    pr3 <- updateCutoffs(pr3, number_of_cutoffs=10)
  
    # check that fdr are the same only when seed is the same
    expect_equal(pr1@fdr, pr2@fdr)
    expect_false(isTRUE(all.equal(pr1@fdr, pr3@fdr)))
})

test_that("updateCutoffs.propd works when ncores > 1", {
    skip(message="WARN:Skipping cannot run on local CI")

    # get propd object and update cutoffs
    set.seed(0)
    pr1 <- propd(x, as.character(y), p=10)
    pr1 <- updateCutoffs(pr1, number_of_cutoffs=10, ncores=1)
    set.seed(0)
    pr2 <- propd(x, as.character(y), p=10)
    pr2 <- updateCutoffs(pr2, number_of_cutoffs=10, ncores=2)

    # check that fdr are the same
    expect_equal(pr1@fdr, pr2@fdr)
})
